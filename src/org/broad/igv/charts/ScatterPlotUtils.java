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

package org.broad.igv.charts;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.collections.DoubleArrayList;

import java.util.*;

/**
 * @author Jim Robinson
 * @date 10/22/11
 */
public class ScatterPlotUtils {

    /**
     * Open a plot on an all loaded data
     * <p/>
     * OTHER, COPY_NUMBER, GENE_EXPRESSION, CHIP, DNA_METHYLATION, TILING_ARRAY, PHASTCON,
     * ALLELE_SPECIFIC_COPY_NUMBER, LOH, MUTATION, RNAI, POOLED_RNAI, CHIP_CHIP, CNV,
     * ALLELE_FREQUENCY, COVERAGE, REPMASK, EXPR
     */
    static HashSet<TrackType> plottableTypes = new HashSet();
    public static final String MUTATION_COUNT = "Mutation Count";

    static {
        plottableTypes.add(TrackType.COPY_NUMBER);
        plottableTypes.add(TrackType.GENE_EXPRESSION);
        plottableTypes.add(TrackType.DNA_METHYLATION);
    }

    public static void openPlot(String chr, int start, int end, int zoom) {

        ScatterPlotData spData = getScatterPlotData(chr, start, end, zoom);
        final ScatterPlotFrame igvPlotFrame = new ScatterPlotFrame(spData);
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                igvPlotFrame.setVisible(true);
            }
        });
    }

    public static boolean hasPlottableTracks() {
        List<Track> tracks = IGV.getInstance().getAllTracks();
        for (Track t : tracks) {
            if (plottableTypes.contains(t.getTrackType())) {
                return true;
            }
        }
        return false;
    }

    private static ScatterPlotData getScatterPlotData(String chr, int start, int end, int zoom) {

        List<Track> tracks = IGV.getInstance().getAllTracks();
        List<String> attributeNames = AttributeManager.getInstance().getAttributeNames();
        LinkedHashMap<String, SampleData> sampleDataMap = new LinkedHashMap<String, SampleData>();
        LinkedHashSet<TrackType> types = new LinkedHashSet<TrackType>();

        // Create a map to keep track of the set of all attribute values for each attribute category.
        LinkedHashMap<String, Set<String>> uniqueAttributeValues = new LinkedHashMap<String, Set<String>>();
        for (String att : attributeNames) {
            uniqueAttributeValues.put(att, new HashSet<String>());
        }
        uniqueAttributeValues.put(MUTATION_COUNT, new HashSet<String>());
        HashSet<String> nonSharedAttributes = new HashSet<String>();

        //String overlayAttribute = IGV.getInstance().getSession().getOverlayAttribute();
        for (Track t : tracks) {

            String sample = t.getSample();
            TrackType type = t.getTrackType();

            if (type == TrackType.MUTATION) {
                // Classify sample by mutation count
                SampleData sampleData = sampleDataMap.get(sample);
                if (sampleData != null) {
                    int mutCount = getMutationCount(chr, start, end, zoom, t);
                    String mutCountString = mutCount < 5 ? String.valueOf(mutCount) : "> 5";
                    sampleData.addAttributeValue(MUTATION_COUNT, mutCountString);
                    uniqueAttributeValues.get(MUTATION_COUNT).add(mutCountString);
                    sampleData.setMutationCount(mutCount);

                }

            } else if (t instanceof DataTrack) {
                DataTrack dataTrack = (DataTrack) t;

                if (plottableTypes.contains(type)) {
                    double regionScore = dataTrack.getAverageScore(chr, start, end, zoom);
                    if (!Double.isNaN(regionScore)) {

                        types.add(type);

                        SampleData sampleData = sampleDataMap.get(sample);
                        if (sampleData == null) {
                            sampleData = new SampleData();
                            sampleDataMap.put(sample, sampleData);
                        }

                        for (String att : attributeNames) {
                            final String attributeValue = dataTrack.getAttributeValue(att);
                            sampleData.addAttributeValue(att, attributeValue);
                            uniqueAttributeValues.get(att).add(attributeValue);
                            final String otherValue = sampleData.getAttributesMap().get(att);
                            if (attributeValue == null) {
                                if (otherValue != null) {
                                    nonSharedAttributes.add(att);
                                }
                            } else if (otherValue == null) {
                                if (attributeValue != null) {
                                    nonSharedAttributes.add(att);
                                }
                            } else {
                                if (!attributeValue.equals(otherValue)) {
                                    nonSharedAttributes.add(att);
                                }
                            }

                        }
                        sampleData.addDataValue(type, regionScore);
                    }
                }

            }
        }

        String[] sampleNames = sampleDataMap.keySet().toArray(new String[sampleDataMap.size()]);

        // Data
        Map<String, double[]> dataMap = new HashMap<String, double[]>(types.size());

        // Loop through track (data) types
        Map<TrackType, Integer> typeCounts = new HashMap<TrackType, Integer>();

        for (TrackType type : types) {
            double[] data = new double[sampleNames.length];
            for (int i = 0; i < sampleNames.length; i++) {
                SampleData sd = sampleDataMap.get(sampleNames[i]);
                // Check for null?  Should be impossible

                double value;
                DoubleArrayList valueList = sd.getData(type);
                if (valueList == null || valueList.isEmpty()) {
                    value = Double.NaN;
                } else if (valueList.size() == 1) {

                    int cnt = typeCounts.get(type) == null ? 1 : typeCounts.get(type) + 1;
                    typeCounts.put(type, cnt);

                    value = valueList.get(0);
                } else {

                    int cnt = typeCounts.get(type) == null ? valueList.size() : typeCounts.get(type) + valueList.size();
                    typeCounts.put(type, cnt);

                    double[] vs = valueList.toArray();
                    value = StatUtils.mean(vs);
                }
                data[i] = value;

            }
            dataMap.put(type.toString(), data);
        }

        // Attributes

        // Get list of "reasonable" attributes with respect to plot series => greater than 1 distinct value, but less
        // than 10.
        List<String> seriesNames = new ArrayList<String>();
        for (Map.Entry<String, Set<String>> entry : uniqueAttributeValues.entrySet()) {
            int cnt = entry.getValue().size();
            String att = entry.getKey();
            if (cnt > 1 && cnt < 10 && !nonSharedAttributes.contains(att) &&
                    !att.equalsIgnoreCase("DATA FILE") && !(att.equalsIgnoreCase("DATA TYPE"))) {
                seriesNames.add(att);
            }
        }

        Map<String, String[]> attMap = new HashMap<String, String[]>(seriesNames.size());

        // Loop through attribute types
        for (String att : seriesNames) {
            String[] attributes = new String[sampleNames.length];
            for (int i = 0; i < sampleNames.length; i++) {
                SampleData sd = sampleDataMap.get(sampleNames[i]);
                // Check for null?  Should be impossible
                final String s = sd.getAttributesMap().get(att);
                attributes[i] = s;
            }
            attMap.put(att, attributes);
        }

        // Collect mutation counts
        int[] mutationCounts = new int[sampleNames.length];
        for (int i = 0; i < sampleNames.length; i++) {
            SampleData sd = sampleDataMap.get(sampleNames[i]);
            // Check for null?  Should be impossible
            mutationCounts[i] = sd.getMutationCount();
        }

        String title = chr + ":" + start + "-" + end;
        return new ScatterPlotData(title, sampleNames, attMap, dataMap, mutationCounts);
    }



    public static int getMutationCount(String chr, int start, int end, int zoom, Track track) {

        return (int) track.getRegionScore(chr, start, end, zoom, RegionScoreType.MUTATION_COUNT, FrameManager.getDefaultFrame().getName());
    }


    /**
     * Container for all data and attributes for a single sample
     */
    public static class SampleData {

        Map<TrackType, DoubleArrayList> valueMap = new HashMap<TrackType, DoubleArrayList>();
        Map<String, String> attributesMap = new HashMap<String, String>();
        int mutationCount;

        public void addDataValue(TrackType type, double value) {

            DoubleArrayList valueArray = valueMap.get(type);
            if (valueArray == null) {
                valueArray = new DoubleArrayList();
                valueMap.put(type, valueArray);
            }
            valueArray.add(value);
        }

        public DoubleArrayList getData(TrackType type) {
            return valueMap.get(type);
        }

        public void addAttributeValue(String att, String attributeValue) {
            attributesMap.put(att, attributeValue);
        }

        public Map<String, String> getAttributesMap() {
            return attributesMap;
        }

        public void setMutationCount(int mutationCount) {
            this.mutationCount = mutationCount;
        }

        public int getMutationCount() {
            return mutationCount;
        }
    }
}
