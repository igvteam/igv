package org.broad.igv.scatterplot;

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

    static {
        plottableTypes.add(TrackType.COPY_NUMBER);
        plottableTypes.add(TrackType.GENE_EXPRESSION);
        plottableTypes.add(TrackType.DNA_METHYLATION);
    }

    public static void openPlot(String chr, int start, int end, int zoom) {

        ScatterPlotData spData = getScatterPlotData(chr, start, end, zoom);
        final org.broad.igv.charts.ScatterPlotFrame igvPlotFrame = new org.broad.igv.charts.ScatterPlotFrame(spData);
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                igvPlotFrame.setVisible(true);
            }
        });

    }

    private static ScatterPlotData getScatterPlotData(String chr, int start, int end, int zoom) {
        List<Track> tracks = IGV.getInstance().getTrackManager().getAllTracks(false);

        List<String> attributeNames = AttributeManager.getInstance().getAttributeNames();
        LinkedHashMap<String, SampleData> sampleDataMap = new LinkedHashMap<String, SampleData>();
        LinkedHashSet<TrackType> types = new LinkedHashSet<TrackType>();
        String key = "samplename";

        LinkedHashMap<String, Set<String>> uniqueAttributeValues = new LinkedHashMap<String, Set<String>>();
        for (String att : attributeNames) {
            uniqueAttributeValues.put(att, new HashSet<String>());
        }
        uniqueAttributeValues.put("Mut Count" , new HashSet<String>());
        HashSet<String> nonSharedAttributes = new HashSet<String>();

        for (Track t : tracks) {
            if (t instanceof DataTrack) {
                DataTrack dataTrack = (DataTrack) t;
                TrackType type = dataTrack.getTrackType();

                if (plottableTypes.contains(type)) {
                    types.add(type);

                    double regionScore = getAverageScore(chr, start, end, zoom, dataTrack);

                    String sample = dataTrack.getAttributeValue(key);

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

            } else if (t.getTrackType() == TrackType.MUTATION) {
                // Classify sample by mutation count
                String sample = t.getAttributeValue(key);
                SampleData sampleData = sampleDataMap.get(sample);
                if (sampleData != null) {
                    int mutCount = getMutationCount(chr, start, end, zoom, t);
                    String mutCountString = mutCount < 5 ? String.valueOf(mutCount) : "> 5";
                    sampleData.addAttributeValue("Mut Count", mutCountString);
                    uniqueAttributeValues.get("Mut Count").add(mutCountString);

                }

            }
        }

        String[] sampleNames = sampleDataMap.keySet().toArray(new String[sampleDataMap.size()]);

        // Data
        Map<String, double[]> dataMap = new HashMap<String, double[]>(types.size());

        // Loop through track (data) types
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
                    value = valueList.get(0);
                } else {
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
            if (cnt > 1 && cnt < 10 && !nonSharedAttributes.contains(att)) {
                seriesNames.add(att);
            }
        }

        Map<String, String[]> attMap = new HashMap<String, String[]>(seriesNames.size());

        for (String att : seriesNames) {
            if (key.equals(att)) continue;
            String[] attributes = new String[sampleNames.length];

            for (int i = 0; i < sampleNames.length; i++) {
                SampleData sd = sampleDataMap.get(sampleNames[i]);
                // Check for null?  Should be impossible

                final String s = sd.getAttributesMap().get(att);
                attributes[i] = s;

            }
            attMap.put(att, attributes);
        }

        return new ScatterPlotData(sampleNames, attMap, dataMap);
    }


    //TODO -- move this to track ?
    private static double getAverageScore(String chr, int start, int end, int zoom, DataTrack dataTrack) {
        double regionScore = 0;
        int intervalSum = 0;
        Collection<LocusScore> scores = dataTrack.getSummaryScores(chr, start, end, zoom);
        for (LocusScore score : scores) {
            if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                int interval = 1; //Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                float value = score.getScore();
                regionScore += value * interval;
                intervalSum += interval;
            }
        }
        if (intervalSum > 0) {
            regionScore /= intervalSum;
        }
        return regionScore;
    }

    private static int getMutationCount(String chr, int start, int end, int zoom, Track track) {

        return (int) track.getRegionScore(chr, start, end, zoom, RegionScoreType.MUTATION_COUNT, FrameManager.getDefaultFrame());
    }


    static class SampleData {

        Map<TrackType, DoubleArrayList> valueMap = new HashMap<TrackType, DoubleArrayList>();
        Map<String, String> attributesMap = new HashMap<String, String>();


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

    }
}
