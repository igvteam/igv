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

package org.broad.igv.tools;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.*;
import org.broad.igv.tools.parsers.CNParser;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.tools.parsers.ToolsWiggleParser;
import org.broad.igv.tools.parsers.UnsortedException;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author jrobinso
 */
public class Preprocessor implements DataConsumer {

    private static Logger log = Logger.getLogger(Preprocessor.class);

    boolean compressed = true;
    private boolean skipZeroes = false;
    private int nZoom = 7;
    int maxExtFactor = 0;
    Zoom[] zoomLevels;
    int nTracks;
    Genome genome;
    Collection<WindowFunction> windowFunctions;
    private String currentChr = "";
    int currentChrLength;
    private int sizeEstimate;
    int nPtsProcessed = 0;
    StatusMonitor statusMonitor;
    double percentComplete = 0.0;
    int lastStartPosition = 0;
    HashSet<String> skippedChromosomes = new HashSet();
    TDFWriter writer;
    Raw rawData;
    Zoom genomeZoom;
    File outputFile;
    ListAccumulator allDataStats;
    List<String> chromosomes = new ArrayList();
    Set<String> visitedChromosomes = new HashSet();
    Map<String, String> attributes = new HashMap();

    List<WindowFunction> allDataFunctions = Arrays.asList(
            WindowFunction.mean,
            WindowFunction.median,
            WindowFunction.min,
            WindowFunction.max,
            WindowFunction.percentile2,
            WindowFunction.percentile10,
            WindowFunction.percentile90,
            WindowFunction.percentile98);


    public Preprocessor(File outputFile,
                        Genome genome,
                        Collection<WindowFunction> windowFunctions,
                        int sizeEstimate,
                        StatusMonitor monitor) {

        this.outputFile = outputFile;
        this.genome = genome;
        this.windowFunctions = windowFunctions;
        this.sizeEstimate = sizeEstimate;
        this.genome = genome;
        this.statusMonitor = monitor;

        allDataStats = new ListAccumulator(allDataFunctions);
        if (statusMonitor == null) {
            PrintStream monStream = this.outputFile == null ? System.err : System.out;
            statusMonitor = new CommandLineStatusMonitor(monStream);
        }
    }

    /**
     * Called to set initial parameters.  It is required that this be called
     * prior to writing the file
     */
    public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames) {
        setTrackParameters(trackType, trackLine, trackNames, true);
    }

    /**
     * Called to set initial parameters.  It is required that this be called
     * prior to writing the file
     */
    public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames, boolean computeWholeGenome) {

        if (trackLine != null) {
            log.info(trackLine);
        }

        if (outputFile != null && writer == null) {
            writer = new TDFWriter(outputFile, genome.getId(), trackType, trackLine, trackNames, windowFunctions, compressed);
            nTracks = trackNames.length;

            // Convert genome coordinates from bp to kbp
            if (computeWholeGenome) {
                int genomeLength = (int) (genome.getNominalLength() / 1000);
                genomeZoom = new Zoom(Globals.CHR_ALL, 0, genomeLength);
            }

            TDFGroup rootGroup = writer.getRootGroup();
            rootGroup.setAttribute("genome", genome.getId());
            rootGroup.setAttribute("maxZoom", String.valueOf(nZoom));

        }
    }


    /**
     * Add an array of data for the given interval.  The array contains a value for each sample/track in this
     * dataset.  The name is an optional probe or feature name.
     */
    public void addData(String chr, int start, int end, float[] data, String name) {

        if (writer == null) {
            return;
        }

        if (skipZeroes) {
            boolean allZeroes = true;
            for (int i = 0; i < data.length; i++) {
                if (data[i] != 0) {
                    allZeroes = false;
                    break;
                }
            }
            if (allZeroes) {
                return;
            }
        }

        // Check for stop signal
        if (statusMonitor != null && statusMonitor.isInterrupted()) {
            throw new PreprocessingException("Preprocessing Halted.");
        }

        if (skippedChromosomes.contains(chr)) {
            return;
        }

        if (currentChr != null && chr.equals(currentChr)) {
            if (start < (lastStartPosition - maxExtFactor)) {
                String msg = "Error: Data is not sorted @ " + chr + " " + start +
                        "  (last position = " + lastStartPosition +
                        "   max ext factor = " + maxExtFactor + ")";
                UnsortedException e = new UnsortedException(msg);
                log.error(msg, e);
                throw e;
            }
        } else {
            newChromosome(chr);
        }

        // Check a second time, in case it just got added
        if (skippedChromosomes.contains(chr)) {
            return;
        }

        // Is this data in range for the chromosome?
        int chrLength = genome.getChromosome(chr).getLength();
        if (start > chrLength) {
            log.warn("Ignoring data from non-existent locus.  Probe = " + name +
                    "  Locus = " + chr + ":" + start + "-" + end + ". " + chr + " length = " + chrLength);
            return;
        }


        // Add to raw data
        rawData.addData(start, end, data, name);

        // Zoom levels
        for (Zoom zl : zoomLevels) {
            zl.addData(start, end, data);
        }

        // Whole genome
        long offset = genome.getCumulativeOffset(chr);
        int gStart = (int) ((offset + start) / 1000);
        int gEnd = Math.max(gStart + 1, (int) ((offset + end) / 1000));


        // Don't include "chrM" in the whole genome view or stats
        if (genomeZoom != null && !(chr.equals("chrM") || chr.equals("M") || chr.equals("MT"))) {
            genomeZoom.addData(gStart, gEnd, data);
            for (int i = 0; i < data.length; i++) {
                allDataStats.add(gEnd - gStart, data[i]);
            }
        }

        lastStartPosition = start;

    }


    /**
     * Start a new chromosome.  Note that data is sorted by chromosome, then start position.
     */
    public void newChromosome(String chr) {

        if (visitedChromosomes.contains(chr)) {
            String msg = "Error: Data is not ordered by start position. Chromosome " + chr +
                    " appears in multiple blocks";
            PreprocessingException e = new PreprocessingException(msg);
            log.error(msg, e);
            throw e;
        }
        visitedChromosomes.add(chr);


        Chromosome c = genome.getChromosome(chr);
        if (c == null) {
            log.warn("Chromosome: " + chr + " not found in .genome file.  Skipping.");
            skippedChromosomes.add(chr);
        } else {

            chromosomes.add(chr);

            log.info("Processing chromosome " + chr);
            if (zoomLevels != null) {
                for (Zoom zl : zoomLevels) {
                    zl.close();
                }
            }
            if (rawData != null) {
                rawData.close();
            }

            currentChr = chr;
            currentChrLength = c.getLength();
            zoomLevels = new Zoom[getNZoom() + 1];
            for (int z = 0; z <= getNZoom(); z++) {
                zoomLevels[z] = new Zoom(chr, z, currentChrLength);
            }

            rawData = new Raw(chr, currentChrLength, 100000);
        }
        lastStartPosition = 0;

    }


    /**
     * Called at end-of-file
     */
    public void parsingComplete() {


    }

    public void finish() {
        if (writer == null) {
            return;
        }

        /**
         * We output the chromosomes in the order specified by the Genome,
         * since that was how the overall genome data was calculated.
         * However we only include those chromosomes which were actually found
         * in the file
         */
        StringBuffer chrString = new StringBuffer();
        Iterator<String> iter = genome.getAllChromosomeNames().iterator();
        while (iter.hasNext()) {
            String chromoName = iter.next();
            if (!chromosomes.contains(chromoName)) {
                continue;
            }
            chrString.append(chromoName);
            if (iter.hasNext()) {
                chrString.append(",");
            }
        }
        writer.getRootGroup().setAttribute(TDFWriter.CHROMOSOMES, chrString.toString());

        for (Map.Entry<String, String> entry : attributes.entrySet()) {
            writer.getRootGroup().setAttribute(entry.getKey(), entry.getValue());
        }

        if (zoomLevels != null) {
            for (Zoom zl : zoomLevels) {
                zl.close();
            }
        }

        if (genomeZoom != null) {
            genomeZoom.close();
        }

        if (rawData == null) {
            // TODO -- delete .tdf file?
            log.warn("No features were found that matched chromosomes in genome: " + genome.getId());
        } else {
            rawData.close();

            // Record max/min
            allDataStats.finish();
            TDFGroup group = writer.getGroup("/");
            group.setAttribute(TDFGroup.USE_PERCENTILE_AUTOSCALING, "true");
            for (WindowFunction wf : allDataFunctions) {
                group.setAttribute(wf.getValue(), String.valueOf(allDataStats.getValue(wf)));
            }
            writer.closeFile();
        }

        if (statusMonitor != null) {
            statusMonitor.setPercentComplete(100);
        }
    }


    public void setType(String type) {
        //this.type = type;
    }

    public void setSortTolerance(int tolerance) {
        maxExtFactor = tolerance;
    }

    public void setAttribute(String key, String value) {
        attributes.put(key, value);
    }

    /**
     * @param sizeEstimate the sizeEstimate to set
     */
    public void setSizeEstimate(int sizeEstimate) {
        this.sizeEstimate = sizeEstimate;
    }

    public void setSkipZeroes(boolean skipZeroes) {
        this.skipZeroes = skipZeroes;
    }

    public int getNZoom() {
        return nZoom;
    }

    public void setNZoom(int nZoom) {
        this.nZoom = nZoom;
    }


    /**
     * Class representing a tile of raw (as opposed to summarized) data.
     */
    class RawTile {
        String dsName;
        int tileNumber;
        int tileStart;
        int tileEnd;
        IntArrayList startArray;
        IntArrayList endArray;
        ArrayList<String> nameList;
        FloatArrayList[] dataArray;

        RawTile(String dsName, int tileNumber, int start, int end) {
            this.dsName = dsName;
            this.tileNumber = tileNumber;
            this.tileStart = start;
            this.tileEnd = end;
            startArray = new IntArrayList();
            endArray = new IntArrayList();
            dataArray = new FloatArrayList[nTracks];
            for (int i = 0; i < nTracks; i++) {
                dataArray[i] = new FloatArrayList();
            }
        }

        void addData(int start, int end, float[] data, String name) {


            if (start > tileEnd) {
                log.warn("Warning: start position > tile end");

            }

            if (end < tileStart) {
                log.warn("Warning: end position > tile end");
            }


            if (name != null && nameList == null) {
                nameList = new ArrayList();
            }


            int dataStart = Math.max(tileStart, start);
            int dataEnd = Math.min(tileEnd, end);
            startArray.add(dataStart);
            endArray.add(dataEnd);
            for (int i = 0; i < data.length; i++) {
                dataArray[i].add(data[i]);
            }
            if (name != null) {
                nameList.add(name);
            }
        }

        void close() {
            try {
                if (startArray.size() > 0) {
                    int[] s = startArray.toArray();
                    int[] e = endArray.toArray();
                    float[][] d = new float[dataArray.length][dataArray[0].size()];
                    for (int i = 0; i < dataArray.length; i++) {
                        d[i] = dataArray[i].toArray();
                    }


                    String[] n = nameList == null ? null : nameList.toArray(new String[]{});
                    TDFBedTile tile = new TDFBedTile(tileStart, s, e, d, n);
                    writer.writeTile(dsName, tileNumber, tile);
                    startArray.clear();
                    endArray.clear();
                    for (int i = 0; i < dataArray.length; i++) {
                        dataArray[i].clear();
                    }
                }
            } catch (IOException ex) {
                log.error(ex.getMessage(), ex);
            }
        }
    }


    /**
     * Class representing the raw dataset
     */
    class Raw {


        String chr;
        String dsName;
        TDFDataset dataset;
        int tileWidth;
        Map<Integer, RawTile> activeTiles = new HashMap();

        Raw(String chr, int chrLength, int tileWidth) {

            this.tileWidth = tileWidth;
            int nTiles = (int) (chrLength / tileWidth) + 1;
            dsName = "/" + chr + "/raw";
            dataset = writer.createDataset(dsName, TDFDataset.DataType.FLOAT, tileWidth, nTiles);

        }

        /**
         * @param start
         * @param end
         * @param data
         */
        public void addData(int start, int end, float[] data, String name) {

            int startTileNumber = (int) (start / tileWidth);
            int endTileNumber = (int) (end / tileWidth);

            // Check for closed tiles -- tiles we are guaranteed to not revisit
            int tmp = (start - maxExtFactor) / tileWidth;
            while (!activeTiles.isEmpty()) {
                Integer tileNumber = activeTiles.keySet().iterator().next();
                if (tileNumber < tmp) {
                    RawTile t = activeTiles.get(tileNumber);
                    t.close();
                    activeTiles.remove(tileNumber);
                } else {
                    break;
                }
            }

            // Add data to all tiles it spans.  The data will be effectively "cut" if it spans multiple tiles.
            for (int t = startTileNumber; t <= endTileNumber; t++) {
                RawTile tile = activeTiles.get(t);
                if (tile == null) {
                    tile = new RawTile(dsName, t, t * tileWidth, (t + 1) * tileWidth);
                    activeTiles.put(t, tile);
                }
                tile.addData(start, end, data, name);
            }

            // Update progress -- assume uniform distribution
            if (statusMonitor != null && sizeEstimate > 0) {
                int p = (int) ((100.0 * nPtsProcessed) / (1.5 * sizeEstimate));
                if (p > percentComplete) {
                    percentComplete = p;
                    statusMonitor.setPercentComplete(percentComplete);
                }
            }
            nPtsProcessed++;
        }

        void close() {
            for (RawTile t : activeTiles.values()) {
                t.close();
            }
            activeTiles = null;
        }


    }

    /**
     * Class representing all the data for a particular zoom level.
     */
    class Zoom {

        int level;
        int tileWidth;
        LinkedHashMap<Integer, Tile> activeTiles = new LinkedHashMap();
        Map<WindowFunction, TDFDataset> datasets = new HashMap();


        Zoom(String chr, int level, int chrLength) {
            int nTiles = (int) Math.pow(2, level);
            tileWidth = chrLength / nTiles + 1;
            this.level = level;

            // Create datasets -- one for each window function
            for (WindowFunction wf : windowFunctions) {
                String dsName = "/" + chr + "/z" + level + "/" + wf.toString();
                datasets.put(wf, writer.createDataset(dsName, TDFDataset.DataType.FLOAT, tileWidth, nTiles));
            }
        }

        public void addData(int start, int end, float[] data) {

            int startTile = start / tileWidth;
            int endTile = end / tileWidth;

            // Check for closed tiles
            int tmp = (start - maxExtFactor) / tileWidth;
            while (!activeTiles.isEmpty()) {
                Integer tileNumber = activeTiles.keySet().iterator().next();
                if (tileNumber < tmp) {
                    Tile t = activeTiles.get(tileNumber);
                    t.close();
                    activeTiles.remove(tileNumber);
                } else {
                    break;
                }
            }


            for (int i = startTile; i <= endTile; i++) {
                Tile t = activeTiles.get(i);
                if (t == null) {
                    t = new Tile(datasets, level, i, 700, tileWidth);
                    activeTiles.put(i, t);
                }
                t.addData(start, end, data);
            }
        }

        // Close all active tiles
        public void close() {
            for (Tile t : activeTiles.values()) {
                t.close();
            }
        }
    }

    /**
     * A tile of summarized data
     */
    class Tile {
        int totalCount;
        int zoomLevel;
        int tileNumber;
        int tileStart;
        int lastFinishedBin = 0;
        double binWidth;
        int nBins;
        int nonEmptyBins;
        ListAccumulator[][] accumulators;
        Map<WindowFunction, TDFDataset> datasets;

        Tile(Map<WindowFunction, TDFDataset> datasets, int zoomLevel, int tileNumber, int nBins, int tileWidth) {
            this.totalCount = 0;
            this.datasets = datasets;
            this.zoomLevel = zoomLevel;
            this.tileNumber = tileNumber;
            this.tileStart = tileNumber * tileWidth;
            this.nBins = nBins;
            this.binWidth = ((double) tileWidth) / nBins;
            this.accumulators = new ListAccumulator[nTracks][nBins];
        }

        /**
         * Add a data point.  Positions are zero based exclusive  (UCSC convention)
         *
         * @param start
         * @param end
         * @param data  array of values at this position,  1 value per track
         */
        void addData(int start, int end, float[] data) {
            totalCount++;

            int startBin = Math.max(0, (int) ((start - tileStart) / binWidth));
            int endBin = Math.min(nBins - 1, (int) ((end - tileStart) / binWidth));

            int tmp = (int) ((start - tileStart - maxExtFactor) / binWidth);


            for (int t = 0; t < nTracks; t++) {
                for (int b = lastFinishedBin; b < tmp; b++) {
                    if (accumulators[t][b] != null) {
                        accumulators[t][b].finish();
                    }
                }
                lastFinishedBin = Math.max(0, tmp - 1);

                for (int b = startBin; b <= endBin; b++) {
                    if (accumulators[t][b] == null) {
                        accumulators[t][b] = new ListAccumulator(datasets.keySet());
                    }
                    accumulators[t][b].add(end - start, data[t]);
                }
            }
        }


        /**
         *
         */
        void close() {

            // Count non-empty bins.  All tracks should be the same
            nonEmptyBins = 0;

            for (int t = 0; t < nTracks; t++) {
                for (int i = 0; i < nBins; i++) {
                    if (accumulators[t][i] != null) {
                        accumulators[t][i].finish();
                        if (t == 0) {
                            nonEmptyBins++;
                        }
                    }
                }
            }

            TDFTile tile = null;

            for (WindowFunction wf : datasets.keySet()) {


                // If < 50% of bins are empty use vary step tile, otherwise use fixed step
                if (nonEmptyBins < 0.5 * nBins) {
                    int[] starts = new int[nonEmptyBins];
                    float[][] data = new float[nTracks][nonEmptyBins];
                    int n = 0;
                    for (int i = 0; i < nBins; i++) {
                        for (int t = 0; t < nTracks; t++) {
                            ListAccumulator acc = accumulators[t][i];
                            if (acc != null) {
                                data[t][n] = acc.getValue(wf);
                                if (t == nTracks - 1) {
                                    starts[n] = (int) (tileStart + (i * binWidth));
                                    n++;
                                }
                            }
                        }
                    }
                    tile = new TDFVaryTile((int) tileStart, binWidth, starts, data);

                } else {
                    float[][] data = new float[nTracks][nBins];
                    for (int t = 0; t < nTracks; t++) {
                        for (int i = 0; i < nBins; i++) {
                            data[t][i] = accumulators[t][i] == null ? Float.NaN : accumulators[t][i].getValue(wf);
                        }
                    }
                    tile = new TDFFixedTile(tileStart, tileStart, binWidth, data);
                }

                String dsName = datasets.get(wf).getName();
                try {
                    writer.writeTile(dsName, tileNumber, tile);
                } catch (IOException exc) {
                    log.error("Error writing tile: " + dsName + " [" + tileNumber + "]", exc);
                    throw new PreprocessingException(exc.getMessage());
                }
            }
        }
    }

    public static boolean isAlignmentFile(String ext) {
        return ext.equalsIgnoreCase(".bam") || ext.equalsIgnoreCase(".sam") ||
                ext.equalsIgnoreCase(".aligned") ||
                ext.equalsIgnoreCase(".bedz") || ext.equalsIgnoreCase(".bed");
    }


    public void preprocess(File iFile, int maxZoomValue, String typeString) throws IOException {

        setNZoom(maxZoomValue);


        String tmp = (typeString == null ? iFile.getAbsolutePath() : typeString).toLowerCase();
        if (tmp.endsWith(".txt")) tmp = tmp.substring(0, tmp.length() - 4);
        if (tmp.endsWith(".gz")) tmp = tmp.substring(0, tmp.length() - 3);
        if (tmp.endsWith("wig") || tmp.endsWith("bedgraph") || tmp.endsWith("cpg") || tmp.endsWith("map")) {
            ToolsWiggleParser wg = new ToolsWiggleParser(iFile.getAbsolutePath(), this, genome);
            wg.parse();
        } else if (tmp.endsWith("cn") || tmp.endsWith("xcn") || tmp.endsWith("igv") || tmp.endsWith("snp")) {
            CNParser cnParser = new CNParser(iFile.getAbsolutePath(), this, genome);
            cnParser.parse();
        } else {
            String msg = "Error: cannot convert files of type '" + tmp + "' to TDF format.";
            msg += "\nTry specifying the file type with the --fileType parameter.";
            throw new RuntimeException(msg);
        }
    }

    /**
     * Strip of "gz" and "txt" extensions to get to true format extension.
     *
     * @param filename
     * @return
     */
    public static String getExtension(String filename) {

        if (filename.endsWith(".gz")) {
            filename = filename.substring(0, filename.length() - 3);
        }

        if (filename.toLowerCase().endsWith(".txt")) {
            filename = filename.substring(0, filename.length() - 4);
        }

        int idx = filename.lastIndexOf('.');
        if (idx < 0) {
            return "";
        } else {
            return filename.substring(idx).toLowerCase();
        }
    }
}
