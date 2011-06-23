/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/**          +-
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Nov 20, 2009
 * Time: 3:43:04 PM
 * To change this template use File | Settings | File Templates.
 */
package org.broad.igv.bbfile;

import org.broad.tribble.util.LittleEndianInputStream;
import org.apache.log4j.Logger;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableFileStream;

import java.io.*;
import java.util.ArrayList;

/*
*   Broad Institute Interactive Genome Viewer Big Binary File (BBFile) Reader
*   -   File reader for UCSC BigWig and BigBed file types.
*
*   Notes:   Table entries refer to Jim Kent of UCSC's document description:
*           "BigWig and BigBed: Enabling Browsing of Large Distributed Data Sets",
*           November 2009.
*
*           The overall binary file layout is defined in Table B of the document.
*
*   BBFile Reader sequences through this binary file layout:
*
*   1) Reads in BBFile Header Table C and determine if file is a valid Big Bed or Big Wig
*      binary file type.
*
*   2) Reads in Zoom Header Tables D if zoom data is present, as defined by zoomLevels
*       in Table C, one for each zoom level.
*
*   3) Reads in  the AutoSQL block Table B if present, as referenced by autoSqlOffset in Table C.
*
*   4) Reads in Total Summary Block Table DD if present, as referenced by
*       TotalSummaryOffset in Table C.
*
*   5) Reads in B+ Tree Header Chromosome Index Header Table E, as referenced
*       by chromosomeTreeOffset in Table C.
*
*   6)Reads in B+ Tree Nodes indexing mChromosome ID's for mChromosome regions;
*       Table F for node type (leaf/child), Table G for leaf items,
*       Table H for child node items.
*
*   7) Reads in R+ Tree Chromosome ID Index Header Table K.
*
*   8) Reads in R+ Tree Nodes indexing of data arranged by mChromosome ID's;
*       Table L for node type (leaf or child), Table M for leaf items,
*       Table N for child node items.
*
*   9) Verifies Data Count of data records, as referenced by fullDataOffset in Table C
*
*   10) References data count records of data size defined in Table M of R+ Tree index
*       for all leaf items in the tree.
*
*   11) Reads in zoom level format Table O for each zoom level comprised of
*       zoom record count followed by that many Table P zoom statistics records,
*       followed by an R+ tree of zoom data locations indexed as in Tables L, M, and N.
*
*   12) Returns information on chromosome name keys and chromosome data regions.
*
*   13) Provides iterators using chromosome names and data regions to extract
*       zoom data, Wig data, and Bed data.
* 
* */

public class BBFileReader {

    public static final long BBFILE_HEADER_OFFSET = 0;

    private static Logger log = Logger.getLogger(BBFileReader.class);

    // Defines the Big Binary File (BBFile) access
    private String path;         // BBFile source file/pathname
    private SeekableStream fis;      // BBFile input stream handle
    private long fileOffset;           // file offset for next item to be read

    private BBFileHeader fileHeader; // Big Binary file header
    private int dataCount;             // Number of data records in the file - Table BB
    private boolean isLowToHigh;       // BBFile binary data format: low to high or high to low
    private int uncompressBufSize;     // buffer byte size for data decompression; 0 for uncompressed

    // AutoSQL String defines custom BigBed formats
    private long autoSqlOffset;
    private String autoSqlFormat;

    // This section defines the zoom items if zoom data exists
    private int zoomLevelCount;       // number of zoom levels defined
    private long zoomLevelOffset;      // file offset to zoom level headers
    private BBZoomLevels zoomLevels;   // zoom level headers and data locations

    // Total Summary Block - file statistical info
    private long mTotalSummaryBlockOffset;
    private BBTotalSummaryBlock totalSummaryBlock;

    // B+ tree
    private long chromIDTreeOffset; // file offset to mChromosome index B+ tree
    private BPTree chromosomeIDTree;     // Container for the mChromosome index B+ tree

    // R+ tree
    private long chromDataTreeOffset;  // file offset to mChromosome data R+ tree
    private RPTree chromosomeDataTree;     // Container for the mChromosome data R+ tree


    public BBFileReader(String path) throws IOException {
        this(path, new SeekableFileStream(new File(path)));

    }

    public BBFileReader(String path, SeekableStream stream) {


        log.debug("Opening BBFile source  " + path);
        this.path = path;
        fis = stream;

        // read in file header
        fileOffset = BBFILE_HEADER_OFFSET;
        fileHeader = new BBFileHeader(path, fis, fileOffset);
        //fileHeader.print();

        if (!fileHeader.isHeaderOK()) {
            log.error("BBFile header is unrecognized type, header magic = " +
                    fileHeader.getMagic());
            throw new RuntimeException("Error reading BBFile header for: " + path);
        }

        // get data characteristics
        isLowToHigh = fileHeader.isLowToHigh();
        uncompressBufSize = fileHeader.getUncompressBuffSize();

        // update file offset past BBFile header
        fileOffset += BBFileHeader.BBFILE_HEADER_SIZE;

        // get zoom level count from file header
        zoomLevelCount = fileHeader.getZoomLevels();

        // extract zoom level headers and zoom data records
        // Note: zoom headers Table D immediately follow the BBFile Header
        if (zoomLevelCount > 0) {

            zoomLevelOffset = fileOffset;

            zoomLevels = new BBZoomLevels(fis, zoomLevelOffset, zoomLevelCount,
                    isLowToHigh, uncompressBufSize);

            // end of zoom level headers - compare with next BBFile item location
            fileOffset += zoomLevelCount * BBZoomLevelHeader.ZOOM_LEVEL_HEADER_SIZE;
        }

        // get the AutoSQL custom BigBed fields
        autoSqlOffset = fileHeader.getAutoSqlOffset();
        if (autoSqlOffset != 0) {
            // read in .as entry
            // mFileOffset = mAutoSqlOffset + sizeof(.as format field);
        }

        // get the Total Summary Block (Table DD)
        fileOffset = fileHeader.getTotalSummaryOffset();
        if (fileHeader.getVersion() >= 2 && fileOffset > 0) {
            totalSummaryBlock = new BBTotalSummaryBlock(fis, fileOffset, isLowToHigh);
            fileOffset += BBTotalSummaryBlock.TOTAL_SUMMARY_BLOCK_SIZE;
        }

        // get Chromosome Data B+ Tree (Table E, F, G, H) : should always exist
        chromIDTreeOffset = fileHeader.getChromosomeTreeOffset();
        if (chromIDTreeOffset != 0) {
            fileOffset = chromIDTreeOffset;
            chromosomeIDTree = new BPTree(fis, fileOffset, isLowToHigh);
        }

        // get number of data records indexed by the R+ chromosome data location tree
        fileOffset = fileHeader.getFullDataOffset();
        dataCount = getDataCount(fis, fileOffset);


    }

    /*
    *   Method returns the Big Binary File pathname.
    *
    *   Returns:
    *       Big Binary File pathname
    * */

    public String getBBFilePath() {
        return path;
    }

    /*
    *   Method returns the Big Binary File input stream handle.
    *
    *   Returns:
    *       Big Binary File input stream handle
    * */

    public SeekableStream getBBFis() {
        return fis;
    }

    /*
    *   Method returns the Big Binary File header which identifies
    *   the file type and content.
    *
    *   Returns:
    *       Big Binary File header (Table C)
    * */

    public BBFileHeader getBBFileHeader() {
        return fileHeader;
    }

    /*
    *   Method returns if the Big Binary File is BigBed.
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is BigBed
    *       (recognized from magic number in file header Table C)
    * */

    public boolean isBigBedFile() {
        return fileHeader.isBigBed();
    }

    /*
    *   Method returns if the Big Binary File is BigWig
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is BigWig
    *       (recognized from magic number in file header Table C)
    * */

    public boolean isBigWigFile() {
        return fileHeader.isBigWig();
    }

    /*
    *   Method returns the total number of data records in the file.
    *
    *   Returns:
    *       Count of the total number of compressed/uncompressed data records:
    *           which for BigBed is the number of bed features,
    *           and for BiGWifg is the number of wig sections.
    * */

    public int getDataCount() {
        return dataCount;
    }

    /*
    *   Method returns the number of chromosomes/contigs in the file.
    *
    *   Note: This is itemCount from B+ tree header BBFile Table E.
    * 
    *   Returns:
    *       Count of the total number of chromosomes/contigs in the file.
    * */

    public long getChromosomeNameCount() {
        return chromosomeIDTree.getItemCount();
    }

    /*
    *   Method returns the number of chromosome/contig regions in the file.
    *
    *   Note: This is itemCount from R+ tree header BBFile Table K.
    *
    *   Returns:
    *       Count of the total number of chromosome/contig regions in the file.
    * */

    public long getChromosomeRegionCount() {
        return chromosomeDataTree.getItemCount();
    }

    /*
   *   Method returns the Big Binary File decompressed buffer size.
   *
   *   Returns:
   *       Largest required buffer size for decompressed data chunks (from Table C)
   * */

    public int getDecompressionBufSize() {
        return uncompressBufSize;
    }

    /*
    *   Method returns if the Big Binary File is written with a low to high byte
    *   order for formatted data.
    *
    *   Returns:
    *       Boolean identifies if Big Binary File is low to high byte order
    *       (recognized from magic number in file header Table C); else is
    *       high to low byte order if false.
    * */

    public boolean isLowToHigh() {
        return isLowToHigh;
    }

    /*
    *   Method returns the total summary block for the Big Binary File.
    *
    *   Returns:
    *       Total summary block data statistics for the whole file (Table DD)
    * */

    public BBTotalSummaryBlock getTotalSummaryBlock() {
        return totalSummaryBlock;
    }

    /*
    *   Method returns the B+ Chromosome Index Tree.
    *
    *   Returns:
    *       B+ Chromosome Index Tree (includes Tables  E, F, G, H)
    * */

    public BPTree getChromosomeIDTree() {
        return chromosomeIDTree;
    }

    /*
    *   Method returns the R+ Chromosome Data Locations Tree.
    *
    *   Returns:
    *       R+ Chromosome Data Locations Tree (includes Tables  K, L, M, N)
    * */

    public RPTree getChromosomeDataTree() {
        return chromosomeDataTree;
    }

    /*
    *   Method returns number of zoom level data is included in the file.
    *
    *   Returns:
    *       Number of zoom levels (from Table C)
    * */

    public int getZoomLevelCount() {
        return zoomLevelCount;
    }

    /*
    *   Method returns the zoom levels in the Big Binary File.
    *
    *   Returns:
    *       Zoom level object containing zoom level headers and R+ zoom data locations tree
    *       (includes Tables  D, O)
    * */

    public BBZoomLevels getZoomLevels() {
        return zoomLevels;
    }

    /*
    *   Method finds the zoom data bounds in R+ tree for a chromosome ID range.
    *
    *   Parameters:
    *       zoomLevel - zoom level
    *       startChromID - start chromosome for the region
    *       endChromID - end chromosome for the region
    *
    *   Returns:
    *       Chromosome region bounds for chromosome ID range
    * */

    public RPChromosomeRegion getZoomLevelBounds(int zoomLevel, int startChromID,
                                                 int endChromID) {

        RPChromosomeRegion chromosomeBounds =
                zoomLevels.getZoomLevelRPTree(zoomLevel).getChromosomeRegion(startChromID, endChromID);

        return chromosomeBounds;
    }

    /*
    *   Method finds chromosome bounds for entire chromosome ID range in the zoom level R+ tree.
    *
    *   Parameters:
    *       zoomLevel - zoom level
    *
    *   Returns:
    *       Chromosome bounds for the entire chromosome ID range in the R+ tree.
    * */

    public RPChromosomeRegion getZoomLevelBounds(int zoomLevel) {

        RPChromosomeRegion chromosomeBounds =
                zoomLevels.getZoomLevelRPTree(zoomLevel).getChromosomeBounds();

        return chromosomeBounds;
    }

    /*
    *   Method returns the zoom record count for the zoom level.
    *
    *   Parameters:
    *       zoomLevel - zoom level
    *
    *   Returns:
    *       Chromosome bounds for the entire chromosome ID range in the R+ tree.
    * */

    public int getZoomLevelRecordCount(int zoomLevel) {

        return zoomLevels.getZoomLevelFormats().get(zoomLevel - 1).getZoomRecordCount();
    }

    /*
    *   Method finds chromosome key name for the associated chromosome ID in the B+ tree.
    *
    *   Returns:
    *       chromosome key name for associated chromosome ID.
    * */

    public String getChromosomeName(int chromID) {

        String chromosomeName = chromosomeIDTree.getChromosomeName(chromID);
        return chromosomeName;
    }

    /*
    *   Method finds chromosome names in the B+ chromosome index tree.
    *
    *   Returns:
    *       LIst of all chromosome key names in the B+ tree.
    * */

    public ArrayList<String> getChromosomeNames() {

        ArrayList<String> chromosomeList = chromosomeIDTree.getChromosomeNames();
        return chromosomeList;
    }

    /*
   *   Returns a chromosome ID  which  can be used to search for a
   *   corresponding data section in the R+ tree for data.
   *
      Parameters:
   *       chromosomeKey - chromosome name of valid key size.
   *
   *
   *   Note: A chromosomeID of -1 means chromosome name not included in B+ tree.
   *
   * */

    public int getChromosomeID(String chromosomeKey) {

        int chromosomeID = chromosomeIDTree.getChromosomeID(chromosomeKey);

        return chromosomeID;
    }

    /*
    *   Method finds the chromosome bounding region in R+ tree for a chromosome ID range.
    *
    *   Parameters:
    *       startChromID - start chromosome for the region
    *       endChromID - end chromosome for the region
    *
    *   Returns:
    *       Chromosome region bounds for chromosome ID range
    * */

    public RPChromosomeRegion getChromosomeBounds(int startChromID, int endChromID) {

        RPChromosomeRegion chromosomeBounds =
                chromosomeDataTree.getChromosomeRegion(startChromID, endChromID);

        return chromosomeBounds;
    }

    /*
    *   Method finds the chromosome bounds for the entire chromosome ID range in the R+ tree.
    *
    *   Returns:
    *       chromosome bounds for the entire chromosome ID range in the R+ tree.
    * */

    public RPChromosomeRegion getChromosomeBounds() {

        RPChromosomeRegion chromosomeBounds = chromosomeDataTree.getChromosomeBounds();

        return chromosomeBounds;
    }

    /*
    *   Method finds all chromosome data regions in the R+ tree.
    *
    *   Returns:
    *       List of chromosome ID's and regions.
    * */

    public ArrayList<RPChromosomeRegion> getChromosomeRegions() {

        ArrayList<RPChromosomeRegion> regionList = chromosomeDataTree.getAllChromosomeRegions();

        return regionList;
    }

    /*
    *   Method finds all zoom level data regions in the R+ tree.
    *
    *   Parameters:
    *       int zoomLevel - zoom level
    *   Returns:
    *       List of chromosome ID's and regions for the zoom level.
    * */

    public ArrayList<RPChromosomeRegion> getZoomLevelRegions(int zoomLevel) {

        ArrayList<RPChromosomeRegion> regionList =
                zoomLevels.getZoomLevelRPTree(zoomLevel).getAllChromosomeRegions();

        return regionList;
    }

    /**
     * Returns an iterator for BigBed features which occupy a chromosome selection region.
     * <p/>
     * Note: the BBFile type should be BigBed; else a null iterator is returned.
     * <p/>
     * Parameters:
     * startChromosome - name of start chromosome
     * startBase     - starting base position for features
     * endChromosome - name of end chromosome
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the specified
     * base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigBed.(see isBigBedFile method)
     */
    public BigBedIterator getBigBedIterator(String startChromosome, int startBase,
                                            String endChromosome, int endBase, boolean contained) {

        if (!isBigBedFile())
            return null;

        if (chromosomeDataTree == null) {
            // get R+ chromosome data location tree (Tables K, L, M, N)
            chromDataTreeOffset = fileHeader.getFullIndexOffset();
            if (chromDataTreeOffset != 0) {
                fileOffset = chromDataTreeOffset;
                chromosomeDataTree = new RPTree(fis, fileOffset, isLowToHigh, uncompressBufSize);
            }

        }
        // go from chromosome names to chromosome ID region
        RPChromosomeRegion selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region
        if (selectionRegion == null)
            throw new RuntimeException("Error finding BigBedIterator region: chromosome not found \n");

        // compose an iterator
        BigBedIterator bedIterator = new BigBedIterator(fis, chromosomeIDTree, chromosomeDataTree,
                selectionRegion, contained);

        return bedIterator;
    }

    /**
     * Returns an iterator for BigBed features for all chromosome regions.
     * <p/>
     * Note: the BBFile type should be BigBed; else a null iterator is returned.
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for all chromosome regions.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigBed.(see isBigBedFile method)
     */
    public BigBedIterator getBigBedIterator() {

        if (!isBigBedFile())
            return null;

        if (chromosomeDataTree == null) {
            // get R+ chromosome data location tree (Tables K, L, M, N)
            chromDataTreeOffset = fileHeader.getFullIndexOffset();
            if (chromDataTreeOffset != 0) {
                fileOffset = chromDataTreeOffset;
                chromosomeDataTree = new RPTree(fis, fileOffset, isLowToHigh, uncompressBufSize);
            }

        }

        // get all region bounds
        RPChromosomeRegion selectionRegion = chromosomeDataTree.getChromosomeBounds();

        // compose an iterator
        boolean contained = true;   /// all regions are contained
        BigBedIterator bedIterator = new BigBedIterator(fis, chromosomeIDTree, chromosomeDataTree,
                selectionRegion, contained);

        return bedIterator;
    }



    /**
     * Returns an iterator for BigWig values which occupy the specified startChromosome region.
     * <p/>
     * Note: the BBFile type should be BigWig; else a null iterator is returned.
     * <p/>
     * Parameters:
     * startChromosome  - name of start chromosome
     * startBase    - starting base position for features
     * endChromosome  - name of end chromosome
     * endBase      - ending base position for feature
     * contained    - flag specifies bed features must be contained in the specified
     * base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigWig.(see isBigWigFile method)
     */
    public BigWigIterator getBigWigIterator(String startChromosome, int startBase,
                                            String endChromosome, int endBase, boolean contained) {

        if (chromosomeDataTree == null) {
            // get R+ chromosome data location tree (Tables K, L, M, N)
            chromDataTreeOffset = fileHeader.getFullIndexOffset();
            if (chromDataTreeOffset != 0) {
                fileOffset = chromDataTreeOffset;
                chromosomeDataTree = new RPTree(fis, fileOffset, isLowToHigh, uncompressBufSize);
            }

        }

        if (!isBigWigFile())
            return null;

        // go from chromosome names to chromosome ID region
        RPChromosomeRegion selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region
        if (selectionRegion == null)
            throw new RuntimeException("Error finding BigWigIterator region: chromosome not found \n");

        // compose an iterator
        BigWigIterator wigIterator = new BigWigIterator(fis, chromosomeIDTree, chromosomeDataTree,
                selectionRegion, contained);

        return wigIterator;
    }

    /**
     * Returns an iterator for BigWig values for all chromosome regions.
     * <p/>
     * Note: the BBFile type should be BigWig; else a null iterator is returned.
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for all chromosome regions.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     * 2) A null object is returned if the file is not BigWig.(see isBigWigFile method)
     */
    public BigWigIterator getBigWigIterator() {

        if (!isBigWigFile())
            return null;

        if (chromosomeDataTree == null) {
            // get R+ chromosome data location tree (Tables K, L, M, N)
            chromDataTreeOffset = fileHeader.getFullIndexOffset();
            if (chromDataTreeOffset != 0) {
                fileOffset = chromDataTreeOffset;
                chromosomeDataTree = new RPTree(fis, fileOffset, isLowToHigh, uncompressBufSize);
            }

        }
        // get all regions bounds
        RPChromosomeRegion selectionRegion = chromosomeDataTree.getChromosomeBounds();

        // compose an iterator
        boolean contained = true;       // all regions are contained
        BigWigIterator wigIterator = new BigWigIterator(fis, chromosomeIDTree, chromosomeDataTree,
                selectionRegion, contained);

        return wigIterator;
    }



    /**
     * Returns an iterator for zoom level records for the chromosome selection region.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * startChromosome - start chromosome name
     * startBase     - staring base position for features
     * endChromosome - end chromosome name
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the
     * specified base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    public ZoomLevelIterator getZoomLevelIterator(int zoomLevel, String startChromosome, int startBase,
                                                  String endChromosome, int endBase, boolean contained) {
        // check for valid zoom level
        if (zoomLevel < 1 || zoomLevel > zoomLevelCount)
            throw new RuntimeException("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree zoomDataTree = zoomLevels.getZoomLevelRPTree(zoomLevel);

        // go from chromosome names to chromosome ID region
        RPChromosomeRegion selectionRegion = getChromosomeBounds(startChromosome, startBase,
                endChromosome, endBase);

        // check for valid selection region  
        if (selectionRegion == null) {
            return ZoomLevelIterator.EmptyIterator.theInstance;
        }

        /// compose an iterator
        ZoomLevelIterator zoomIterator = new ZoomLevelIterator(fis, chromosomeIDTree,
                zoomDataTree, zoomLevel, selectionRegion, contained);

        return zoomIterator;
    }

    /**
     * Returns an iterator for zoom level records for all chromosome regions.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    public ZoomLevelIterator getZoomLevelIterator(int zoomLevel) {

        // check for valid zoom level
        if (zoomLevel < 1 || zoomLevel > zoomLevelCount)
            throw new RuntimeException("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree zoomDataTree = zoomLevels.getZoomLevelRPTree(zoomLevel);

        // get all regions bounds
        RPChromosomeRegion selectionRegion = zoomDataTree.getChromosomeBounds();

        // compose an iterator
        boolean contained = true;   //all regions are contained
        ZoomLevelIterator zoomIterator = new ZoomLevelIterator(fis, chromosomeIDTree,
                zoomDataTree, zoomLevel, selectionRegion, contained);

        return zoomIterator;
    }

    /**
     * Returns an iterator for zoom level records for the chromosome selection region.
     * <p/>
     * Note: the BBFile can be BigBed or BigWig.
     * <p/>
     * Parameters:
     * zoomLevel - zoom level for data extraction; levels start at 1
     * selectionRegion - chromosome selection region consists of:
     * startChromID - ID of starting chromosome
     * startBase     - staring base position for features
     * endChromID - ID of endind chromosome
     * endBase       - ending base position for feature
     * contained     - flag specifies bed features must be contained in the
     * specified base region if true; else can intersect the region if false
     * <p/>
     * Returns:
     * Iterator to provide BedFeature(s) for the requested chromosome region.
     * Error conditions:
     * 1) An empty iterator is returned if region has no data available
     */
    public ZoomLevelIterator getZoomLevelIterator(int zoomLevel, RPChromosomeRegion selectionRegion,
                                                  boolean contained) {
        // check for valid zoom level
        if (zoomLevel < 1 || zoomLevel > zoomLevelCount)
            throw new RuntimeException("Error: ZoomLevelIterator zoom level is out of range\n");

        // get the appropriate zoom level R+ zoom data index tree
        RPTree zoomDataTree = zoomLevels.getZoomLevelRPTree(zoomLevel);

        /// compose an iterator
        ZoomLevelIterator zoomIterator = new ZoomLevelIterator(fis, chromosomeIDTree,
                zoomDataTree, zoomLevel, selectionRegion, contained);

        return zoomIterator;
    }

    /*
    *   Method generates a chromosome bounds region for the supplied chromosome region name.
    *
    *   Note: No attempt is made to verify the region exists in the file data, nor
    *   which data is being examined.
    *
    *   Parameters:
    *       startChromosome - name of start chromosome
    *       startBase - starting base position for region
    *       endChromosome - name of end chromosome
    *       endBase - ending base position for region
    *
    *   Returns:
    *       Chromosome bounds of a named chromosome region for data extraction;
    *       or null for regions not found in the B+ chromosome index tree.
    * */

    private RPChromosomeRegion getChromosomeBounds(String startChromosome, int startBase,
                                                   String endChromosome, int endBase) {

        // find the chromosome ID's using the name to get a valid name key, then associated ID
        String startChromKey = chromosomeIDTree.getChromosomeKey(startChromosome);
        int startChromID = chromosomeIDTree.getChromosomeID(startChromKey);
        if (startChromID < 0)       // mChromosome not in data?
            return null;

        String endChromKey = chromosomeIDTree.getChromosomeKey(endChromosome);
        int endChromID = chromosomeIDTree.getChromosomeID(endChromKey);
        if (endChromID < 0)       // mChromosome not in data?
            return null;

        // create the bounding mChromosome region
        RPChromosomeRegion chromBounds = new RPChromosomeRegion(startChromID, startBase,
                endChromID, endBase);

        return chromBounds;
    }

    /*
    *   Method reads data count which heads the data section of the BBFile.
    *
    *   Returns:
    *       Data count of the number of data records:
    *          number of Bed features for BigBed
    *          number of Wig sections for BigWig
    * */

    private int getDataCount(SeekableStream fis, long fileOffset) {
        int dataCount;
        LittleEndianInputStream lbdis = null;
        DataInputStream bdis = null;

        // Note: dataCount in BBFile is simply a 4 byte int
        // positioned at fullDataOffset in Table C
        byte[] buffer = new byte[4];

        try {
            // read dataCount into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);

            // decode data count with proper byte stream reader
            // first assume byte order is low to high
            if (isLowToHigh) {
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
                dataCount = lbdis.readInt();
            } else {
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));
                dataCount = bdis.readInt();
            }
        } catch (IOException ex) {
            throw new RuntimeException("Error reading data count for all data", ex);
        }

        // data count was read properly
        return dataCount;
    }

} // end of BBFileReader
