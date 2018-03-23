/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.LRUCache;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * @author jrobinso
 */
public class TDFReader {

    static final Logger log = Logger.getLogger(TDFReader.class);
    public static final int GZIP_FLAG = 0x1;

    private SeekableStream seekableStream = null;
    private int version;
    private Map<String, IndexEntry> datasetIndex;
    private Map<String, IndexEntry> groupIndex;
    private TrackType trackType;
    private String trackLine;
    private String[] trackNames;
    private String genomeId;
    LRUCache<String, TDFGroup> groupCache = new LRUCache(20);
    LRUCache<String, TDFDataset> datasetCache = new LRUCache(20);
    TDFTile wgTile;

    Map<WindowFunction, Double> valueCache = new HashMap();
    private List<WindowFunction> windowFunctions;
    ResourceLocator locator;

    boolean compressed = false;

    Set<String> chrNames;
    private final CompressionUtils compressionUtils;

    //private String path;

    public static TDFReader getReader(String path) {
        return getReader(new ResourceLocator(path));
    }

    public static TDFReader getReader(ResourceLocator locator) {
        return new TDFReader(locator);
    }

    public TDFReader(ResourceLocator locator) {
        //this.path = path;
        this.locator = locator;
        try {
            log.debug("Getting stream");
            seekableStream = IGVSeekableStreamFactory.getInstance().getStreamFor(locator.getPath());
            log.debug("Reading header");
            readHeader();
            log.debug("Done reading header");

        } catch (IOException ex) {
            log.error("Error loading file: " + locator.getPath(), ex);
            throw new DataLoadException("Error loading file: " + ex.toString(), locator.getPath());
        }
        compressionUtils = new CompressionUtils();
    }

    public void close() {
        try {
            seekableStream.close();
        } catch (IOException e) {
            log.error("Error closing reader for: " + getPath(), e);
        }
    }

    public String getPath() {
        return locator.getPath();
    }

    private void readHeader() throws IOException {

        // Buffer for the magic number, version, index position, and index
        // byte count + header byte count  (4 + 4 + 8 + 4 + 4)
        //byte[] buffer = new byte[24];
        //readFully(buffer);
        byte[] buffer = readBytes(0, 24);
        ByteBuffer byteBuffer = ByteBuffer.wrap(buffer);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        int magicNumber = byteBuffer.getInt();

        byte[] magicBytes = new byte[4];
        System.arraycopy(buffer, 0, magicBytes, 0, 4);
        String magicString = new String(magicBytes);

        if (!(magicString.startsWith("TDF") || magicString.startsWith("IBF"))) {
            String msg = "Error reading header: bad magic number.";
            throw new DataLoadException(msg, locator.getPath());
        }

        version = byteBuffer.getInt();
        long idxPosition = byteBuffer.getLong();
        int idxByteCount = byteBuffer.getInt();
        int nHeaderBytes = byteBuffer.getInt();

        byte[] bytes = readBytes(24, nHeaderBytes);
        byteBuffer = ByteBuffer.wrap(bytes);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        if (version >= 2) {
            int nWFs = byteBuffer.getInt();
            log.debug("nWFs: " + nWFs);
            this.windowFunctions = new ArrayList(nWFs);
            for (int i = 0; i < nWFs; i++) {
                String wfName = StringUtils.readString(byteBuffer);
                try {
                    windowFunctions.add(WindowFunction.valueOf(wfName));
                } catch (Exception e) {
                    log.error("Error creating window function: " + wfName, e);
                }
            }
        } else {
            windowFunctions = Arrays.asList(WindowFunction.mean);
        }

        // Track type
        try {
            trackType = TrackType.valueOf(StringUtils.readString(byteBuffer));

        } catch (Exception e) {
            trackType = TrackType.OTHER;
        }

        // Track line
        trackLine = StringUtils.readString(byteBuffer).trim();


        int nTracks = byteBuffer.getInt();
        trackNames = new String[nTracks];
        for (int i = 0; i < nTracks; i++) {
            trackNames[i] = StringUtils.readString(byteBuffer);
        }

        // Flags
        if (version > 2) {
            genomeId = StringUtils.readString(byteBuffer);
            int flags = byteBuffer.getInt();
            compressed = (flags & GZIP_FLAG) != 0;
        } else {
            compressed = false;
        }


        readMasterIndex(idxPosition, idxByteCount);

    }


    private void readMasterIndex(long idxPosition, int nBytes) throws IOException {

        try {
//fis.seek(idxPosition);
            //byte[] bytes = new byte[nBytes];
            //readFully(bytes);
            byte[] bytes = readBytes(idxPosition, nBytes);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

            int nDatasets = byteBuffer.getInt();

            datasetIndex = new LinkedHashMap(nDatasets);
            for (int i = 0; i < nDatasets; i++) {
                String name = StringUtils.readString(byteBuffer);
                long fPosition = byteBuffer.getLong();
                int n = byteBuffer.getInt();
                datasetIndex.put(name, new IndexEntry(fPosition, n));
            }

            int nGroups = byteBuffer.getInt();
            groupIndex = new LinkedHashMap(nGroups);
            for (int i = 0; i < nGroups; i++) {
                String name = StringUtils.readString(byteBuffer);
                long fPosition = byteBuffer.getLong();
                int n = byteBuffer.getInt();
                groupIndex.put(name, new IndexEntry(fPosition, n));
            }
        } catch (BufferUnderflowException e) {
            // We intermittently see this exception in this method.  Log as much info as possible
            log.error("BufferUnderflowException.  path=" + getPath() + "  idxPosition=" + idxPosition + "  nBytes=" + nBytes);
            throw e;
        }
    }

    public TDFDataset getDataset(String chr, int zoom, WindowFunction windowFunction) {

        // Version 1 only had mean
        String wf = getVersion() < 2 ? "" : "/" + windowFunction.toString();
        String zoomString = chr.equals(Globals.CHR_ALL) ? "0" : String.valueOf(zoom);
        String dsName = "/" + chr + "/z" + zoomString + wf;

        TDFDataset ds = getDataset(dsName);


        return ds;
    }

    public synchronized TDFDataset getDataset(String name) {

        if (datasetCache.containsKey(name)) {
            return datasetCache.get(name);
        }

        try {
            if (datasetIndex.containsKey(name)) {
                IndexEntry ie = datasetIndex.get(name);
                long position = ie.position;
                int nBytes = ie.nBytes;

                //fis.seek(position);
                //byte[] buffer = new byte[nBytes];
                //readFully(buffer);
                byte[] buffer = readBytes(position, nBytes);
                ByteBuffer byteBuffer = ByteBuffer.wrap(buffer);
                byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

                TDFDataset ds = new TDFDataset(name, byteBuffer, this);
                datasetCache.put(name, ds);
                return ds;
            } else {
                datasetCache.put(name, null);
                return null;
            }

        } catch (IOException ex) {
            log.error("Error reading dataset: " + getPath() + " (" + name + ")", ex);
            throw new RuntimeException("System error occured while reading dataset: " + name);
        }
    }

    public Collection<String> getDatasetNames() {
        return datasetIndex.keySet();
    }

    public Collection<String> getGroupNames() {
        return groupIndex.keySet();
    }

    public synchronized TDFGroup getGroup(String name) {
        if (groupCache.containsKey(name)) {
            return groupCache.get(name);
        }

        try {
            IndexEntry ie = groupIndex.get(name);
            long position = ie.position;
            int nBytes = ie.nBytes;

            //fis.seek(position);
            //byte[] buffer = new byte[nBytes];
            //readFully(buffer);
            byte[] buffer = readBytes(position, nBytes);
            ByteBuffer byteBuffer = ByteBuffer.wrap(buffer);
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

            TDFGroup group = new TDFGroup(name, byteBuffer);

            groupCache.put(name, group);

            return group;

        } catch (IOException ex) {
            log.error("Error reading group: " + name, ex);
            throw new RuntimeException("System error occured while reading group: " + name);
        }
    }

    // TODO -- move to dataset class

    public TDFTile readTile(TDFDataset ds, int tileNumber) {

        try {
            if (tileNumber >= ds.tilePositions.length) {
                // TODO - return empty tile
                return null;
            }

            long position = ds.tilePositions[tileNumber];
            if (position < 0) {

                // Indicates empty tile
                // TODO -- return an empty tile?
                return null;
            }

            int nBytes = ds.tileSizes[tileNumber];
            //fis.seek(position);
            //byte[] buffer = new byte[nBytes];
            //readFully(buffer);
            byte[] buffer = readBytes(position, nBytes);
            if (compressed) {
                buffer = compressionUtils.decompress(buffer);

            }

            return TileFactory.createTile(buffer, trackNames.length);
        } catch (IOException ex) {
            String tileName = ds.getName() + "[" + tileNumber + "]";
            log.error("Error reading data tile: " + tileName, ex);
            throw new RuntimeException("System error occured while reading tile: " + tileName);
        }
    }

    /**
     * @return the version
     */
    public int getVersion() {
        return version;
    }

    /**
     * @return the trackType
     */
    public TrackType getTrackType() {
        return trackType;
    }

    /**
     * @return the trackLine
     */
    public String getTrackLine() {
        return trackLine;
    }

    /**
     * @return the trackNames
     */
    public String[] getTrackNames() {
        return trackNames;
    }

    /**
     * Only available for version >= 3.
     *
     * @return
     */
    public String getGenomeId() {
        return genomeId;
    }

    private Double getValue(WindowFunction wf) {
        if (!valueCache.containsKey(wf)) {
            TDFGroup rootGroup = getGroup("/");
            String maxString = rootGroup.getAttribute(wf.getValue());
            try {
                valueCache.put(wf, Double.parseDouble(maxString));
            } catch (Exception e) {
                log.info("Warning: value '" + wf.toString() + "' not found in tdf value " + getPath());
                valueCache.put(wf, null);
            }
        }
        return valueCache.get(wf);
    }

    /**
     * Return the default upper value for the data range.  A check is made to see if both upper and lower limits
     * are equal to zero, within a tolerance.  If they are the upper limit is arbitrarily set to "100".  This is
     * protection against the pathological case that can occur with datasets consisting of largely zero values.
     *
     * @return
     */
    public double getUpperLimit() {
        Double val = getValue(WindowFunction.percentile98);
        double upperLimit = val == null ? getDataMax() : val;
        if (upperLimit < 1.0e-30 && getLowerLimit() < 1.0e-30) {
            upperLimit = 100;
        }
        return upperLimit;
    }

    public double getLowerLimit() {
        Double val = getValue(WindowFunction.percentile2);
        return val == null ? getDataMin() : val;
    }

    public double getDataMax() {
        Double val = getValue(WindowFunction.max);
        return val == null ? 100 : val;
    }

    public double getDataMin() {
        Double val = getValue(WindowFunction.min);
        return val == null ? 0 : val;
    }


    public synchronized byte[] readBytes(long position, int nBytes) throws IOException {
        seekableStream.seek(position);
        byte[] buffer = new byte[nBytes];
        int read = seekableStream.read(buffer, 0, nBytes);
        return buffer;
    }

    /**
     * @return the windowFunctions
     */
    public List<WindowFunction> getWindowFunctions() {
        return windowFunctions;
    }

    /**
     * Return a list of all chromosome names represented in this dataset
     * TODO -- record this explicitly in the TDF file
     *
     * @return
     */

    public Set<String> getChromosomeNames() {
        if (chrNames == null) {
            ///DatasetIndex /chr1/z0/mean=org.broad.igv.tdf.TDFReader$IndexEntry@6a493b65
            chrNames = new HashSet();
            for (String key : datasetIndex.keySet()) {
                String[] tokens = Globals.forwardSlashPattern.split(key);
                int nTokens = tokens.length;
                if (nTokens > 1) {
                    chrNames.add(tokens[1]);
                }
            }
        }
        return chrNames;
    }

    class IndexEntry {

        long position;
        int nBytes;

        public IndexEntry(long position, int nBytes) {
            this.position = position;
            this.nBytes = nBytes;
        }
    }


    /**
     * Print index entries (for debugging)
     */
    public void dumpIndex() {

        for (Map.Entry<String, IndexEntry> entry : datasetIndex.entrySet()) {

            String dsName = entry.getKey();

            TDFDataset ds = getDataset(dsName);

            int size = 0;
            for (int sz : ds.tileSizes) {
                size += sz;
            }
            System.out.println(dsName + "\t" + size);

            datasetCache.clear();
        }
    }

    public TDFTile getWholeGenomeTile(Genome genome, WindowFunction wf) {

        if (wgTile == null) {

            int binCount = 700;
            int nTracks = this.getTrackNames().length; // TODO -- is there a more direct way to know this?
            double binSize = (genome.getNominalLength() / 1000) / binCount;
            Accumulator[][] accumulators = new Accumulator[nTracks][binCount];

            for (String chrName : genome.getLongChromosomeNames()) {

                TDFDataset chrDataset = getDataset(chrName, 0, wf);
                if(chrDataset == null) continue;

                List<TDFTile> chrTiles = chrDataset.getTiles();
                chrDataset.clearCache(); // Don't cache these
                for (TDFTile t : chrTiles) {
                    int[] chrStart = t.getStart();
                    int[] chrEnd = t.getEnd();

                    for (int p = 0; p < chrStart.length; p++) {

                        int gStart = genome.getGenomeCoordinate(chrName, chrStart[p]);
                        int gEnd = genome.getGenomeCoordinate(chrName, chrEnd[p]);

                        int binStart = (int) (gStart / binSize);
                        int binEnd = Math.min(binCount - 1, (int) (gEnd / binSize));
                        for (int b = binStart; b <= binEnd; b++) {

                            for (int n = 0; n < nTracks; n++) {

                                Accumulator acc = accumulators[n][b];
                                if (acc == null) {
                                    acc = new Accumulator(WindowFunction.mean);
                                    accumulators[n][b] = acc;
                                }

                                int basesCovered = Math.min(gEnd, binEnd) - Math.max(gStart, binStart);
                                acc.add(basesCovered, t.getData(n)[p], null);
                            }
                        }
                    }
                }
            }

            // Fill data array for "wg tile"
            float[][] data = new float[nTracks][binCount];
            for (int n = 0; n < nTracks; n++) {
                Accumulator[] accArray = accumulators[n];
                for (int p = 0; p < accArray.length; p++) {
                    Accumulator acc = accArray[p];
                    if (acc != null) {
                        acc.finish();
                        data[n][p] = acc.getValue();

                    }
                }
            }

            //    public TDFFixedTile(int tileStart, int start, double span, float[][] data) {
            wgTile = new TDFFixedTile(0, 0, binSize, data);

        }
        return wgTile;
    }
}