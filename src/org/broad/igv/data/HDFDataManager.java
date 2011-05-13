/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
/*
 * FeatureTileManager.java
 *
 * Created on November 13, 2007, 9:57 PM
 *
 * To change this template, choose Tools | Template Manager
 * and openFile the template in the editor.
 */
package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.h5.DataAccessException;
import org.broad.igv.h5.HDF5Reader;
import org.broad.igv.h5.HDF5ReaderFactory;
import org.broad.igv.h5.ObjectNotFoundException;
import org.broad.igv.session.RendererFactory;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.*;

import java.util.*;

/**
 * Manages data from a single HDF5 file.
 *
 * @author jrobinso
 */
public class HDFDataManager {

    private static Logger log = Logger.getLogger(HDFDataManager.class);
    private static Integer rawIndex = -1;
    private ObjectCache<String, SummaryTile2D> summaryTileCache = new ObjectCache();
    private ObjectCache<String, int[]> intArrayCache = new ObjectCache();
    private Map<String, float[]> dataMaxCache = new HashMap(25);
    private Map<String, float[]> dataMinCache = new HashMap(25);
    private Map<String, Double> rawTileSizeCache = new HashMap();
    private Map<String, int[]> rawIndecesCache = new HashMap();
    private Map<String, Integer> longestFeatureCache = new HashMap();
    private Map<String, Integer> numZoomLevels = new HashMap();
    private Map<String, Integer> chrLengths = new HashMap();
    private Map<String, Double> binWidthCache = new HashMap();
    /**
     * Cache of zoom level parameters.  Key is the zoom level
     */
    private LRUCache<String, ZoomParameters> zoomParameterCache = new LRUCache(this, 100);
    /**
     * Cache of rawMeansPath arrays -- shared across all instances.  Used for location arrays
     */
    private ObjectCache<String, Object> locationsCache = new ObjectCache(100);
    private HDF5Reader dataReader;
    /**
     * The track type.  Read from HDF5 file
     */
    TrackType trackType = TrackType.OTHER;
    /**
     * The data window size in base pairs.  Used for fix width  data
     */
    int windowSpan = 1;
    /**
     * Unit (in bp) for location values,  so 1000 == kb,  10000000 = mb, etc
     */
    int locationUnit = 1;
    /**
     * Flag indicating if this file has a data section
     */
    private boolean hasData = true;
    /**
     * Flag indication if feature end locations are included in the file
     */
    private boolean hasEnd = true;
    /**
     * Path prefix for data
     */
    private String dataGroupPath;
    /**
     * Path prefix for location
     */
    private String featureGroupPath;
    private String[] trackNames;
    private boolean normalized;
    private int version;
    private WindowFunction windowFunction;
    private TrackProperties trackProperties;
    private ResourceLocator locator;
    private Genome genome;

    /**
     * Creates a new instance of FeatureTileManager for a given track and chromosome
     *
     * @param locator
     * @throws DataAccessException
     */

    public HDFDataManager(ResourceLocator locator) throws DataAccessException {
        this(locator, null);
    }

    public HDFDataManager(ResourceLocator locator, Genome genome) throws DataAccessException {
        this.genome = genome;
        int rootGroup = -1;
        try {
            setGroupPaths();

            this.locator = locator;
            dataReader = HDF5ReaderFactory.getReader(locator);

            // TODO -- open file and read these
            // this.h5File = dataReader.openFile(locator.getPath());
            rootGroup = dataReader.openGroup("/");
            readRootAttributes(rootGroup);
            readTrackProperties(rootGroup);

            this.windowFunction = isCopyNumber() ? WindowFunction.median : WindowFunction.mean;

            // TODO read "type" attribute
            // boolean hasData = dataReader.readIntegerAttribute(rootGroup, "has.data") != 0;

        } finally {
            if (rootGroup > 0) {
                dataReader.closeGroup(rootGroup);
            }
        }

    }

    public boolean isCopyNumber() {
        TrackType tt = getTrackType();
        return tt == TrackType.COPY_NUMBER || tt == TrackType.ALLELE_SPECIFIC_COPY_NUMBER;
    }


    public ResourceLocator getResourceLocator() {
        return locator;
    }

    /**
     * Method description
     *
     * @return
     */
    public String[] getTrackNames() {
        if (trackNames == null) {
            int datasetId = -1;
            try {
                String dsName = "/data/track.id";
                datasetId = dataReader.openDataset(dsName);

                // A hack for ChIP data.
                List<String> names = dataReader.readAllStrings(datasetId);
                trackNames = new String[names.size()];
                if (this.getTrackType() == TrackType.CHIP) {
                    String strip = ".aligned";
                    for (int i = 0; i < names.size(); i++) {
                        String nm = names.get(i);
                        trackNames[i] = nm.endsWith(strip)
                                ? nm.substring(0, nm.length() - strip.length()) : nm;
                    }
                } else {
                    for (int i = 0; i < names.size(); i++) {
                        trackNames[i] = names.get(i);
                    }
                }

            } finally {
                if (datasetId > 0) {
                    dataReader.closeDataset(datasetId);
                }
            }
        }
        return trackNames;
    }

    /**
     * Return the number of precomputed zoom levels for a chromosome
     *
     * @param chr
     * @return
     */
    int getNumZoomLevels(String chr) {
        Integer nz = numZoomLevels.get(chr);
        if (nz == null) {

            if (chr.equals("All")) {
                nz = 1;
            } else {
                int featureGroup = dataReader.openGroup(featureGroupPath + chr);
                nz = dataReader.readIntegerAttribute(featureGroup, "zoom.levels");
                dataReader.closeGroup(featureGroup);
            }

            numZoomLevels.put(chr, nz);
        }

        return nz.intValue();

    }

    /**
     * Return the length of a chromosome.
     * TODO -- not sure why this is here
     *
     * @param chr
     * @return
     */
    public int getChrLength(String chr) {
        Integer len = chrLengths.get(chr);
        if (len == null) {
            if (chr.equals("All")) {
                len = 1;
            } else {
                int featureGroup = dataReader.openGroup(featureGroupPath + chr);
                len = dataReader.readIntegerAttribute(featureGroup, "length");
                dataReader.closeGroup(featureGroup);
            }

            chrLengths.put(chr, len);
        }

        return len.intValue();
    }

    /**
     * @return the track type for this rawMeansPath
     */
    public TrackType getTrackType() {
        return trackType;
    }

    /**
     * Read the summarized data values from the rawMeansPath for all tracks over
     * the given locus range.
     *
     * @param dsName
     * @param startIndex
     * @param endIndex
     * @return
     */
    float[][] readData(String dsName, int startIndex, int endIndex) {
        int datasetId = -1;
        try {
            datasetId = dataReader.openDataset(dsName);
            float[][] dataArray = dataReader.readDataSlice(datasetId, startIndex, endIndex);
            return dataArray;
        } finally {
            if (datasetId > 0) {
                dataReader.closeDataset(datasetId);
            }

        }
    }

    /**
     * Method description
     *
     * @param trackNumber
     * @param chr
     * @return
     */
    public double getDataMax(int trackNumber, String chr) {

        // Try track properties first
        if (trackProperties != null && !Float.isNaN(trackProperties.getMaxValue())) {
            return trackProperties.getMaxValue();
        }

        String rawMeansPath = dataGroupPath + chr + "/raw/mean";
        String rawStdPath = dataGroupPath + chr + "/raw/stddev";
        String key = rawMeansPath;
        float[] tmp = dataMaxCache.get(key);
        if (tmp == null) {

            int meanDatasetId = -1;
            int stdDatasetId = -1;
            try {
                meanDatasetId = dataReader.openDataset(rawMeansPath);
                float[] means = dataReader.readAllFloats(meanDatasetId);
                stdDatasetId = dataReader.openDataset(rawStdPath);
                float[] stds = dataReader.readAllFloats(stdDatasetId);

                if (means.length == 0 || stds.length == 0) {
                    return 10;
                }
                tmp = new float[means.length];
                for (int i = 0; i < means.length; i++) {
                    tmp[i] = means[i] + 2 * stds[i];
                }
            } catch (Exception objectNotFoundException) {
                return 10.0f;
            } finally {
                if (meanDatasetId > 0) {
                    dataReader.closeDataset(meanDatasetId);
                }
                if (stdDatasetId > 0) {
                    dataReader.closeDataset(stdDatasetId);
                }

            }
        }
        return tmp[trackNumber];

    }

    /**
     * Method description
     *
     * @param trackNumber
     * @param chr
     * @return
     */
    public double getDataMin(int trackNumber, String chr) {

        // Try track properties first
        if (trackProperties != null && !Float.isNaN(trackProperties.getMinValue())) {
            return trackProperties.getMinValue();
        }

        String dataset = dataGroupPath + chr + "/raw/min";
        String key = dataset;
        float[] tmp = dataMinCache.get(key);
        if (tmp == null) {
            int datasetId = -1;
            try {
                datasetId = dataReader.openDataset(dataset);
                tmp = dataReader.readAllFloats(datasetId);
                dataMinCache.put(key, tmp);
            } catch (ObjectNotFoundException objectNotFoundException) {
                return 0.0f;
            } finally {
                if (datasetId > 0) {
                    dataReader.closeDataset(datasetId);
                }

            }
        }

        if (tmp.length == 0) {
            return 0;
        }

        // TODO -- hack to get around bogus min values recorded by chip-seq processor
        double min = tmp[trackNumber];
        if (min < 30000) {
            min = 0;
        }
        return min;
    }

    protected double getAttributeForGroup(String groupName, String attribute) {
        int group = dataReader.openGroup(groupName);
        double att = dataReader.readDoubleAttribute(group, attribute);

        dataReader.closeGroup(group);

        return att;
    }

    protected void setGroupPaths() {
        dataGroupPath = "/data/";
        featureGroupPath = "/features/";
    }

    /**
     * Return the zoomm dependent parameters for the given level.
     */
    private ZoomParameters getZoomParameters(int zoom, String chr) {
        String zoomGroupPath = featureGroupPath + chr + ((zoom == rawIndex) ? "/raw" : "/z" + zoom);
        String key = zoomGroupPath;
        ZoomParameters params = zoomParameterCache.get(key);

        if (params == null) {
            int zoomGroup = dataReader.openGroup(zoomGroupPath);
            double tileWidth = dataReader.readDoubleAttribute(zoomGroup, "tile.width");
            String h5Name = zoomGroupPath + "/tile.boundary";

            int datasetId = dataReader.openDataset(h5Name);
            int[] tileBoundaries = dataReader.readAllInts(datasetId);
            dataReader.closeDataset(datasetId);

            params = new ZoomParameters();
            params.zoom = zoom;
            params.tileWidth = tileWidth;
            params.tileBoundaries = tileBoundaries;

            if (zoom != rawIndex) {
                double maxCount = dataReader.readDoubleAttribute(zoomGroup, "max.count");
                double meanCount = dataReader.readDoubleAttribute(zoomGroup, "mean.count");
                double medianCount = dataReader.readDoubleAttribute(zoomGroup, "median.count");

                params.maxCount = maxCount;
                params.medianCount = medianCount;
                params.meanCount = meanCount;
                synchronized (zoomParameterCache) {
                    zoomParameterCache.put(key, params);
                }
            }

            dataReader.closeGroup(zoomGroup);
        }

        return params;
    }

    /**
     * Estimate the mean value over all occupied bins of counts per bin.
     * If the zoom level is precomputed (zoom < zoomLevels) this value is
     * recorded in the HDF file and is exact.  If the zoom level is not it.
     * We can't know the exact value since the number of occupied bins is
     * not known.
     *
     * @param zoom
     * @return
     */
    private double estimateMeanCount(int zoom, String chr) {
        int nz = getNumZoomLevels(chr);

        if (zoom < nz) {
            return getZoomParameters(zoom, chr).meanCount;
        } else {
            ZoomParameters zp = getZoomParameters(nz - 1, chr);

            if (zp != null) {
                double factor = Math.pow(2, zoom - nz + 1);

                return Math.max(1.0, zp.meanCount / factor);
            } else {

                // TODO -- handle file without any precomputed zoom levels ?
                return 1.0;
            }

        }
    }

    protected List<SummaryTile2D> getPrecomputedSummaryTiles(String chr, int startLocation,
                                                             int endLocation, int zoom, WindowFunction windowFunction) {

        ZoomParameters zParams = getZoomParameters(zoom, chr);
        int startTile = (int) (startLocation / zParams.tileWidth);
        int endTile = (int) Math.min(zParams.tileBoundaries.length - 1,
                endLocation / zParams.tileWidth);
        int nTiles = endTile - startTile + 1;
        if (nTiles <= 0) {
            log.debug(
                    "Negative tile count.  ntiles = " + nTiles + " startLoction = " + startLocation + " endLocation = " + endLocation + " startTile = " + startTile + " endTile = " + endTile);
            return new ArrayList();
        }

        List<SummaryTile2D> tiles = new ArrayList(nTiles);
        for (int t = 0; t < nTiles; t++) {
            int tileNumber = startTile + t;
            int tileStart = (int) (tileNumber * zParams.tileWidth);
            SummaryTile2D summaryTile = getSummaryTile(chr, tileNumber, tileStart, zoom);
            if (summaryTile != null) {
                tiles.add(summaryTile);
            }

        }
        return tiles;
    }

    protected SummaryTile2D getSummaryTile(String chr, int tileNumber, int startLocation, int zoom) {
        String key = "" + chr + "_" + zoom + "_" + tileNumber + "_" + windowFunction;
        SummaryTile2D summaryTile = summaryTileCache.get(key);

        if (summaryTile == null) {
            SummaryDataTile tile = loadSummaryDataTile(chr, tileNumber, zoom);

            if (tile.isEmpty()) {
                return null;
            }

            double meanCount = estimateMeanCount(zoom, chr);

            summaryTile = new SummaryTile2D(tileNumber, startLocation);
            for (int i = 0; i < tile.getSize(); i++) {
                for (int j = 0; j < tile.getNTracks(); j++) {
                    float value = tile.getValue(j, i);

                    if (!Float.isNaN(value)) {

                        summaryTile.addScore(j, new SummaryScore(
                                (tile.getStart(i) * locationUnit),
                                (int) (tile.getEnd(i) * locationUnit), value));
                    }

                }
            }
            summaryTileCache.put(key, summaryTile);
        }

        return summaryTile;
    }

    private SummaryDataTile loadSummaryDataTile(String chr, int tileNumber, int zoom) {
        ZoomParameters zParams = getZoomParameters(zoom, chr);
        int startIndex = (tileNumber == 0) ? 0 : zParams.tileBoundaries[tileNumber - 1];
        int endIndex = zParams.tileBoundaries[tileNumber];

        // Get the bin width
        double binWidth = 1;
        try {
            String key = featureGroupPath + chr + "/z" + zoom;
            Double tmp = binWidthCache.get(key);
            if (tmp == null) {
                int featureGroup = dataReader.openGroup(key);
                tmp = new Double(dataReader.readDoubleAttribute(featureGroup, "bin.size"));
                binWidthCache.put(key, tmp);
                dataReader.closeGroup(featureGroup);
            }
            binWidth = tmp.floatValue();
        } catch (ObjectNotFoundException e) {
        }


        SummaryDataTile tile = new SummaryDataTile(tileNumber, binWidth);

        // Check for an empty tile
        if (startIndex == endIndex) {
            return tile;
        }

        if (startIndex > endIndex) {
            System.err.println("Start index > end index!");

            return tile;
        }


        try {
            float[][] data = null;

            if ((windowFunction == null) || (windowFunction == WindowFunction.count)) {
                String dsName = featureGroupPath + chr + "/z" + zoom + "/" + "count";
                int datasetId = dataReader.openDataset(dsName);
                float[] counts = dataReader.readFloats(datasetId, startIndex, endIndex);
                dataReader.closeDataset(datasetId);
                tile.setCounts(counts);

            }
            if ((windowFunction != null) && (windowFunction != WindowFunction.count)) {
                String dsName = dataGroupPath + chr + "/z" + zoom + "/" + windowFunction;
                try {
                    data = readData(dsName, startIndex, endIndex);

                } catch (Exception exception) {
                    exception.printStackTrace();
                }
            }

            String dsName = featureGroupPath + chr + "/z" + zoom + "/" + "start";
            String cacheKey = dsName + "_" + tileNumber;
            int[] startLocations = (int[]) locationsCache.get(cacheKey);
            if (startLocations == null) {
                int ds = dataReader.openDataset(dsName);
                startLocations = dataReader.readInts(ds, startIndex, endIndex);
                dataReader.closeDataset(ds);
                locationsCache.put(cacheKey, startLocations);
            }

            tile.setLocations(startLocations);
            tile.setValues(data);

        } catch (Exception ex) {
            log.error(
                    "Warning: problem loading data tile: z= " + zoom + " start= " + startIndex + " end= " + endIndex,
                    ex);
        }

        return tile;
    }

    /**
     * Method description
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @return
     */
    public List<SummaryTile2D> getSummaryTilesForRange(String chr, int startLocation,
                                                       int endLocation, int zoom) {

        // assert endLocation >= startLocation;

        int startTile = 0;
        int endTile = 0;
        if (zoom < getNumZoomLevels(chr)) {
            return getPrecomputedSummaryTiles(chr, startLocation, endLocation, zoom,
                    windowFunction);
        } else {
            int chrLength = getChrLength(chr);
            endLocation = Math.min(endLocation, chrLength);

            // Get some extra tiles to be sure we get long features that start
            // before startLocation
            int longestFeature = getLongestFeature(chr);

            int adjustedStart = Math.max(0, startLocation - longestFeature);
            int adjustedEnd = Math.min(chrLength, endLocation);

            int tmp = Math.max(15, this.getNumZoomLevels(chr) + 1);
            int z = Math.min(tmp, zoom);
            int nTiles = (int) Math.pow(2, z);
            double tileWidth = ((double) chrLength) / nTiles;


            startTile = (int) (adjustedStart / tileWidth);
            endTile = (int) (Math.min(chrLength, adjustedEnd) / tileWidth);

            List<SummaryTile2D> tiles = new ArrayList(nTiles);
            for (int tileNumber = startTile; tileNumber <= endTile; tileNumber++) {
                int tileStart = (int) (tileNumber * tileWidth);
                int tileEnd = Math.min(chrLength, (int) ((tileNumber + 1) * tileWidth));

                String key = chr + "_" + zoom + "_" + tileNumber + "_" + windowFunction;
                SummaryTile2D tile = summaryTileCache.get(key);
                if (tile == null) {
                    tile = computeSummaryTile(chr, tileNumber, tileStart, tileEnd);
                    if (tile != null) {
                        synchronized (summaryTileCache) {
                            summaryTileCache.put(key, tile);
                        }
                    }
                }
                tiles.add(tile);
            }

            return tiles;
        }

    }

    /**
     * Return the longest feature in the rawMeansPath for the given chromosome.  This
     * is needed when computing summary data for a region.
     * <p/>
     * TODO - This default implementaiton is crude and should be overriden by subclasses.
     *
     * @param chr
     * @return
     */
    private synchronized int getLongestFeature(String chr) {

        if (!longestFeatureCache.containsKey(chr)) {
            int maxLongestFeature = 1000;

            if (genome != null) {
                String genomeId = genome.getId();
                // Put a limit on the "longest feature" to prevent outliers from dictating
                // huge data loads.  This is neccessary due to a flaw in the indexing scheme.
                maxLongestFeature = 5000;
            }
            int longestFeature = 1000;

            String groupName = featureGroupPath + chr + "/raw";
            int groupId = dataReader.openGroup(groupName);
            try {
                longestFeature = dataReader.readIntegerAttribute(groupId, "longest.feature");
            } catch (Exception e) {
                if (getTrackType() != TrackType.GENE_EXPRESSION && !hasEnd) {
                    longestFeature = this.windowSpan;
                }

            }
            longestFeature = Math.min(longestFeature, maxLongestFeature);
            longestFeatureCache.put(chr, longestFeature);

        }
        return longestFeatureCache.get(chr);
    }

    /**
     * Compute a tile of summary scores from raw data.
     * TODO -- this method is very similar to the like named method in AbstractDataSource.  Rationalize
     * and combine the two.
     *
     * @param chr
     * @param tileNumber
     * @param startLocation
     * @param endLocation
     * @return
     */
    private SummaryTile2D computeSummaryTile(String chr, int tileNumber, int startLocation,
                                             int endLocation) {

        // The binWidth is the width in base pairs of a single screen pixel
        // TODO -- binWidth should be passed in
        double binWidth = FrameManager.getDefaultFrame().getScale();

        // Get the raw data as a list of tiles.
        List<DataTile2D> tiles = getRawData(chr, startLocation, endLocation);

        // Create a list of features (scores) to process, one per sample (i.e. values row).
        List<List<LocusScore>> smallFeatures = new ArrayList();

        // Loop through the raw data.  If the data interval is > the bin
        // width (i.e > a screen pixel)  use it as is.  If it is < a screen pixel
        // put it in the list of "smallFeatures" for combining with other scores
        // tat fall on this pixel.
        SummaryTile2D tile = new SummaryTile2D(tileNumber, startLocation);
        for (DataTile2D rawTile : tiles) {
            if ((rawTile != null) && !rawTile.isEmpty()) {

                int[] starts = rawTile.getStartLocations();
                int[] ends = rawTile.getEndLocations();
                float[][] values = rawTile.getValues();
                List<List<LocusScore>> scoresToSegregate = new ArrayList();
                List<List<LocusScore>> scoresToBin = new ArrayList();

                for (int i = 0; i < values.length; i++) {
                    scoresToSegregate.add(new ArrayList(1000));
                    scoresToBin.add(new ArrayList(1000));
                }

                for (int i = 0; i < starts.length; i++) {
                    int start = starts[i];
                    int end = (ends == null) ? start + 1 : Math.max(start + 1, ends[i]);
                    for (int trackNumber = 0; trackNumber < values.length; trackNumber++) {
                        float val = values[trackNumber][i];

                        // NaN is treated as a missing value (as if its not there)
                        if (!Float.isNaN(val)) {
                            SummaryScore ss = new SummaryScore(start, end, val);
                            double pixelWidth = (end - start) / binWidth;
                            if (pixelWidth < 1) {
                                scoresToBin.get(trackNumber).add(ss);
                            } else {
                                scoresToSegregate.get(trackNumber).add(ss);
                            }
                        }

                    }
                }

                for (int trackNumber = 0; trackNumber < scoresToSegregate.size(); trackNumber++) {
                    List<LocusScore> locusScores = scoresToSegregate.get(trackNumber);

                    // Disabled - for now, segregatng overlaping features
                    // ProcessingUtils.segregateScores(scoresToSegregate.get(trackNumber),
                    // windowFunction);


                    for (LocusScore score : locusScores) {
                        tile.addScore(trackNumber, score);
                    }

                    List<FeatureBin> bins = computeBins(startLocation, endLocation,
                            scoresToBin.get(trackNumber));
                    for (FeatureBin fBin : bins) {

                        int binStart = fBin.getStart();
                        int binEnd = (int) (binStart + binWidth);
                        int counts = fBin.getFeatureCount();
                        float binValue = 0;

                        if ((windowFunction == null) || (windowFunction == WindowFunction.count)) {
                            binValue = counts;
                        } else {
                            float[] tmp = fBin.getFeatureScores();

                            if (tmp == null) {
                                binValue = Float.NaN;
                            } else {
                                binValue = ProcessingUtils.computeStat(tmp, windowFunction);
                            }

                        }

                        // For now the interpretation of "Float.NaN" is missing data.  No reason
                        // to add these.
                        if (!Float.isNaN(binValue)) {
                            tile.addScore(trackNumber,
                                    new SummaryScore(binStart, binEnd, binValue));
                        }
                    }

                }


            }
        }

        if (tile != null) {
            for (List<LocusScore> scores : tile.getSummaryScores().values()) {

                if (scores != null) {
                    FeatureUtils.sortFeatureList(scores);
                }

            }
        }
        return tile;
    }

    private static List<FeatureBin> computeBins(int startLocation, int endLocation,
                                                List<LocusScore> features) {
        double binSize = FrameManager.getDefaultFrame().getScale();
        int nBins = (int) Math.ceil((endLocation - startLocation) / binSize);
        double correctedBinSize = ((double) (endLocation - startLocation)) / nBins;
        return (new FeatureBinCalculator()).computeFeatureBins(features, nBins, correctedBinSize,
                startLocation, endLocation);
    }

    /**
     * @param trackNumber
     * @param tiles
     * @return
     * @deprecated
     */
    private List<LocusScore> aggregateScores(int trackNumber, List<SummaryTile2D> tiles) {
        List<LocusScore> joinedScores = new ArrayList(tiles.size() * 700);
        LocusScore previousScore = null;

        for (SummaryTile2D tile : tiles) {
            for (LocusScore score : tile.getScores(trackNumber)) {
                if (!Float.isNaN(score.getScore())) {

                    if (previousScore == null) {
                        previousScore = new SummaryScore(score);
                        joinedScores.add(previousScore);
                    } else {
                        if (score.getScore() == previousScore.getScore()) {

                            // score values are identical, stretch the current score
                            // The only purpose of this is to reconsitute segmented data from
                            // cn files that were created from cbs (segmented) files.
                            previousScore.setEnd(score.getEnd());

                        } else {

                            SummaryScore newScore = new SummaryScore(score);

                            // score value has changed. Adjust end of previous
                            // score (if any), and start of this score to meet 1/2
                            // way
                            int delta = newScore.getStart() - previousScore.getEnd();
                            previousScore.setEnd(previousScore.getEnd() + delta / 2);
                            newScore.setStart(previousScore.getEnd());

                            joinedScores.add(newScore);

                            previousScore = newScore;
                        }

                    }
                }

            }
        }


        return joinedScores;
    }    /*
     *    private void dumpFeatures(List<IGVFeature> features) {
     * PrintWriter pw = null;
     * try {
     * pw = new PrintWriter(new FileWriter("feature_dump.txt"));
     * for (IGVFeature f : features) {
     * pw.println(f.getStart());
     * }
     * } catch (IOException ex) {
     * Logger.getLogger(HDFDataManager.class.getName()).log(Level.SEVERE, null, ex);
     * } finally {
     * pw.closeFile();
     * }
     * }
     */


    /**
     * Method description
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    public List<IGVFeature> getFeaturesForRange(String chr, int startLocation, int endLocation) {
        return null;
    }

    ObjectCache<String, DataTile2D> rawDataCache = new ObjectCache(100);

    /**
     * Retrieve a list of raw data tiles covering the specified locus range.
     * <p/>
     * This method is package scope (as opposed to private) to permit unit testing.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    List<DataTile2D> getRawData(String chr, int startLocation, int endLocation) {


        List<DataTile2D> tiles = new ArrayList();

        // Tiles are defined to be a constant size in base pairs.  Get the size
        // frm the hdf5 file.
        if (!rawTileSizeCache.containsKey(chr)) {
            String groupName = featureGroupPath + chr + "/raw";
            int groupId = dataReader.openGroup(groupName);
            synchronized (rawTileSizeCache) {
                rawTileSizeCache.put(chr, dataReader.readDoubleAttribute(groupId, "index.span"));
            }

        }

        double rawTileSize = rawTileSizeCache.get(chr);

        // Get the start indeces array.  This has the index of the first data point
        // for each tile.   Since tiles are constant width in base pairs they will
        // contain a varying number of data points.
        if (!rawIndecesCache.containsKey(chr)) {
            String dsName = featureGroupPath + chr + "/raw/index";
            synchronized (rawIndecesCache) {
                rawIndecesCache.put(chr, readAllInts(dsName));
            }
        }

        int[] rawIndeces = rawIndecesCache.get(chr);

        // Calculate the start tile number
        int startTile = Math.max(0, (int) (startLocation / rawTileSize) - 1);
        if (startTile >= rawIndeces.length) {
            List<DataTile2D> tmp = new ArrayList();
            tmp.add(DataTile2D.getNullDataTile());
            return tmp;
        }

        int startIndex = rawIndeces[startTile];

        // Calculate the end tile number
        int endTile = (int) (endLocation / rawTileSize) + 1;
        int endIndex = (endTile >= rawIndeces.length)
                ? rawIndeces[rawIndeces.length - 1] : rawIndeces[endTile];

        // Equal start and end indeces indicates that no data falls in this range.
        if (endIndex <= startIndex) {
            tiles.add(DataTile2D.getNullDataTile());
            return tiles;
        } else {

            String keyPrefix = "raw_" + chr + "_";
            while (startTile < endTile) {
                String key = keyPrefix + startTile;
                DataTile2D tile = rawDataCache.get(key);
                if (tile == null) {
                    break;
                }

                tiles.add(tile);
                startTile++;

            }

            // If there are any tiles not found in the cache load the remainder.
            // In theory we could check for tiles at the end, but if a data load
            // is required getting some extra data will not hurt.
            if (startTile == endTile) {
                return tiles;
            }

            // Update start index to that of the first non-empty tile
            startIndex = rawIndeces[startTile];

            // If startIndes == endIndex the remaining tiles have no data
            if (endIndex <= startIndex) {
                for (int t = startTile; t < endTile; t++) {
                    tiles.add(DataTile2D.getNullDataTile());
                }
            } else {

                String dsName = featureGroupPath + chr + "/raw/start";
                int[] starts = readInts(dsName, startIndex, endIndex);

                // end is optional
                int[] ends = null;
                if (hasEnd) {
                    try {
                        dsName = featureGroupPath + chr + "/raw/end";
                        ends = readInts(dsName, startIndex, endIndex);
                    } catch (Exception ex) {

                        // This is expected in normal operation.  Nothing to do.
                        hasEnd = false;
                    }

                }

                if (ends == null) {
                    ends = new int[starts.length];
                    for (int i = 0; i < starts.length; i++) {
                        ends[i] = starts[i] + windowSpan;
                    }

                }


                // Try to get data
                // TODO hasData should be recorded and read from the file as an
                // attrbute.
                float[][] data = null;

                if (hasData) {
                    dsName = dataGroupPath + chr + "/raw/value";
                    try {
                        data = this.readData(dsName, startIndex, endIndex);    // dataReader.readFloats(h5File, dsName, startIndex, endIndex);
                    } catch (ObjectNotFoundException ex) {
                        log.info("No data object found for: " + dsName);
                        hasData = false;
                        List<DataTile2D> tmp = new ArrayList();
                        tmp.add(DataTile2D.getNullDataTile());
                        return tmp;

                    }

                }

                List<DataTile2D> newTiles = createDataTiles(startTile, endTile, rawIndeces, starts,
                        ends, data);
                tiles.addAll(newTiles);

                int tileNumber = startTile;
                for (DataTile2D tile : newTiles) {
                    String key = keyPrefix + tileNumber;
                    synchronized (rawDataCache) {
                        rawDataCache.put(key, tile);
                    }
                    tileNumber++;

                }
            }


        }
        return tiles;
    }

    /**
     * Method description
     *
     * @param startTile
     * @param endTile
     * @param indeces
     * @param startLocations
     * @param endLocations
     * @param data
     * @return
     */
    public List<DataTile2D> createDataTiles(int startTile, int endTile, int[] indeces,
                                            int[] startLocations, int[] endLocations, float[][] data) {

        List<DataTile2D> tiles = new ArrayList();
        if (indeces.length <= startTile) {
            return tiles;
        }

        int nTracks = data.length;
        int i0 = indeces[startTile];

        for (int tile = startTile; tile < endTile; tile++) {
            int startIdx = ((tile >= indeces.length) ? startLocations.length - 1 : indeces[tile]);
            int endIdx = ((tile + 1 >= indeces.length)
                    ? startLocations.length - 1 : indeces[tile + 1]);
            if (endIdx <= startIdx) {
                tiles.add(DataTile2D.getNullDataTile());
            } else {

                int nPts = endIdx - startIdx;
                int[] s = new int[nPts];
                int[] e = new int[nPts];
                float[][] d = new float[nTracks][nPts];

                int idx = startIdx - i0;

                System.arraycopy(startLocations, idx, s, 0, nPts);
                System.arraycopy(endLocations, idx, e, 0, nPts);
                for (int t = 0; t < nTracks; t++) {
                    System.arraycopy(data[t], idx, d[t], 0, nPts);
                }

                tiles.add(new DataTile2D(s, e, d));

            }

        }
        return tiles;

    }

    int[] readInts(String dsName, int startIndex, int endIndex) {
        String key = dsName + "_" + startIndex + "_" + endIndex;
        int[] dataArray = intArrayCache.get(key);
        if (dataArray == null) {
            int datasetId = dataReader.openDataset(dsName);
            dataArray = dataReader.readInts(datasetId, startIndex, endIndex);
            dataReader.closeDataset(datasetId);
            intArrayCache.put(key, dataArray);
        }

        return dataArray;
    }

    int[] readAllInts(String dsName) {
        String key = dsName + "_ALL";
        int[] dataArray = intArrayCache.get(key);
        if (dataArray == null) {
            int datasetId = dataReader.openDataset(dsName);
            dataArray = dataReader.readAllInts(datasetId);
            dataReader.closeDataset(datasetId);
            intArrayCache.put(key, dataArray);
        }

        return dataArray;
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isNormalized() {
        return normalized;
    }

    private void readRootAttributes(int rootGroup) throws DataAccessException {


        try {
            version = dataReader.readIntegerAttribute(rootGroup, "version");
        } catch (ObjectNotFoundException objectNotFoundException) {
            throw new DataAccessException("Unsupported HDF5 version number: 0");
        }
        if (version < 2) {
            throw new DataAccessException("Unsupported HDF5 version number: " + version);
        }

        try {
            normalized = (dataReader.readIntegerAttribute(rootGroup, "normalized")) > 0;
        } catch (Exception e) {

            // Default
            normalized = false;
        }
        try {
            trackType = TrackType.valueOf(dataReader.readStringAttribute(rootGroup,
                    "type").toUpperCase());
        } catch (Exception e) {

            // Default
            trackType = TrackType.OTHER;
        }

        try {
            windowSpan = dataReader.readIntegerAttribute(rootGroup, "window.span");
        } catch (ObjectNotFoundException e) {
            windowSpan = 1;
        }


        try {
            locationUnit = dataReader.readIntegerAttribute(rootGroup, "locationUnit");
        } catch (Exception e) {

            // Default
            locationUnit = 1;
        }

        // TODO read "type" attribute
        // boolean hasData = dataReader.readIntegerAttribute(rootGroup, "has.data") != 0;


    }

    private void readTrackProperties(int group) {

        trackProperties = new TrackProperties();
        boolean foundProperties = false;
        String propertyString = null;
        try {
            propertyString = dataReader.readStringAttribute(group, "trackline");
        } catch (ObjectNotFoundException e) {
            // This is normal, nothing to do.
        }
        if (propertyString != null && propertyString.trim().length() > 10) {
            foundProperties = ParsingUtils.parseTrackLine(propertyString, trackProperties);
        }

        if (!foundProperties) {
            // Old method for reading properties
            try {
                String altColorString = dataReader.readStringAttribute(group, "track.altColor");
                getTrackProperties().setAltColor(
                        ColorUtilities.stringToColor(altColorString));
            } catch (Exception e) {
            }
            try {
                String midColorString = dataReader.readStringAttribute(group, "track.midColor");
                getTrackProperties().setMidColor(
                        ColorUtilities.stringToColor(midColorString));
            } catch (Exception e) {
            }
            try {
                String colorString = dataReader.readStringAttribute(group, "track.color");
                getTrackProperties().setColor(ColorUtilities.stringToColor(colorString));
            } catch (Exception e) {
            }
            try {
                String description = dataReader.readStringAttribute(group, "track.description");
                getTrackProperties().setDescription(description);
            } catch (Exception e) {
            }
            try {
                int height = dataReader.readIntegerAttribute(group, "track.height");
                getTrackProperties().setHeight(height);
            } catch (Exception e) {
            }
            try {
                String autoscale = dataReader.readStringAttribute(group, "track.autoscale");
                getTrackProperties().setAutoScale(Boolean.valueOf(autoscale));
            } catch (Exception e) {
            }
            try {
                float minValue = (float) dataReader.readDoubleAttribute(group, "track.minValue");
                getTrackProperties().setMinValue(minValue);
            } catch (Exception e) {
            }
            try {
                float midValue = (float) dataReader.readDoubleAttribute(group, "track.midValue");
                getTrackProperties().setMidValue(midValue);
            } catch (Exception e) {
            }
            try {
                float maxValue = (float) dataReader.readDoubleAttribute(group, "track.maxValue");
                getTrackProperties().setMaxValue(maxValue);
            } catch (Exception e) {
            }
            try {
                String drawMidValue = dataReader.readStringAttribute(group, "track.drawMidValue");
                getTrackProperties().setDrawYLine(Boolean.valueOf(drawMidValue));
            } catch (Exception e) {
            }
            try {
                String name = dataReader.readStringAttribute(group, "track.name");
                getTrackProperties().setName(name);
            } catch (Exception e) {
            }
            try {
                String rendererClass = dataReader.readStringAttribute(group, "track.renderer");
                getTrackProperties().setRendererClass(RendererFactory.getRendererClass(rendererClass));
            } catch (Exception e) {
            }
            try {
                String windowFunction = dataReader.readStringAttribute(group, "track.windowFunction");
                getTrackProperties().setWindowingFunction(WindowFunction.valueOf(windowFunction));
            } catch (Exception e) {
            }
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    public void setWindowFunction(WindowFunction windowFunction) {
        this.windowFunction = windowFunction;
    }

    class SummaryDataTile {

        /**
         * The tileNumber.  Used to support unit tests, and  general  debugging
         */
        private int tileNumber;

        // window span in bp
        private int binWidth = 1;
        int[] locations;
        private float[] counts;
        float[][] values;

        /**
         * Constructs ...
         *
         * @param tileNumber
         * @param binWidth
         */
        public SummaryDataTile(int tileNumber, double binWidth) {
            this.tileNumber = tileNumber;
            this.binWidth = (int) Math.round(binWidth);
        }

        /**
         * Method description
         *
         * @return
         */
        public int getNTracks() {
            return (values == null) ? 0 : values.length;
        }

        /**
         * Method description
         *
         * @param locations
         */
        public void setLocations(int[] locations) {
            this.locations = locations;
        }

        /**
         * Method description
         *
         * @param scores
         */
        public void setValues(float[][] scores) {
            this.values = scores;
        }

        /**
         * Method description
         *
         * @return
         */
        public int getSize() {
            assert values.length == locations.length;
            return locations.length;
        }

        /**
         * Method description
         *
         * @return
         */
        public boolean isEmpty() {
            return (locations == null) || (locations.length == 0);
        }

        /**
         * @param index
         * @return
         */
        public int getStart(int index) {
            return locations[index];
        }

        /**
         * @param index
         * @return
         */
        public float getEnd(int index) {
            return locations[index] + binWidth;
        }

        /**
         * @param trackNumber
         * @param index
         * @return
         */
        public float getValue(int trackNumber, int index) {
            return (values == null) ? Float.NaN : values[trackNumber][index];
        }

        /**
         * Method description
         *
         * @param index
         * @return
         */
        public float getCount(int index) {
            return counts[index];
        }

        /**
         * Method description
         *
         * @return
         */
        public int getTileNumber() {
            return tileNumber;
        }

        /**
         * Method description
         *
         * @return
         */
        public float[] getCounts() {
            return counts;
        }

        /**
         * Method description
         *
         * @param counts
         */
        public void setCounts(float[] counts) {
            this.counts = counts;
        }
    }

    class ZoomParameters {

        private int zoom;
        private int[] tileBoundaries;
        private double tileWidth;
        private double maxCount;
        private double meanCount;
        private double medianCount;
    }
}
