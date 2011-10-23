package org.broad.igv.scatterplot;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
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

        List<Track> tracks = IGV.getInstance().getTrackManager().getAllTracks(false);

        List<String> attributeNames = AttributeManager.getInstance().getAttributeNames();
        LinkedHashMap<String, SampleData> sampleDataMap = new LinkedHashMap<String, SampleData>();
        LinkedHashSet<TrackType> types = new LinkedHashSet<TrackType>();

        for (Track t : tracks) {
            if (t instanceof DataTrack) {
                DataTrack dataTrack = (DataTrack) t;
                TrackType type = dataTrack.getTrackType();
                     if(type == TrackType.GENE_EXPRESSION) {
                    System.out.println();
                }
                if (plottableTypes.contains(type)) {
                    types.add(type);

                    double regionScore = getAverageScore(chr, start, end, zoom, dataTrack);

                    String key = "LINKING_ID";
                    String sample = dataTrack.getAttributeValue(key);

                    SampleData sampleData = sampleDataMap.get(sample);
                    if (sampleData == null) {
                        sampleData = new SampleData();
                        for (String att : attributeNames) {
                            sampleData.addAttributeValue(att, dataTrack.getAttributeValue(att));
                        }
                        sampleDataMap.put(sample, sampleData);
                    }
                    sampleData.addDataValue(type, regionScore);


                }
            }
        }

        String[] sampleNames = sampleDataMap.keySet().toArray(new String[sampleDataMap.size()]);

        // Data
        Map<String, double[]> dataMap = new HashMap<String, double[]>(types.size());
        for (TrackType type : types) {
            double[] data = new double[sampleNames.length];
            for (int i = 0; i < sampleNames.length; i++) {
                SampleData sd = sampleDataMap.get(sampleNames[i]);
                // Check for null?  Should be impossible

                double value;
                DoubleArrayList valueList = sd.getData(type);
                if (valueList == null || valueList.isEmpty()) {
                    value = Double.NaN;
                }
                else if (valueList.size() == 1) {
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
        Map<String, String[]> attMap = new HashMap<String, String[]>(attributeNames.size());
        for (String att : attributeNames) {
            String[] attributes = new String[sampleNames.length];
            for (int i = 0; i < sampleNames.length; i++) {
                SampleData sd = sampleDataMap.get(sampleNames[i]);
                // Check for null?  Should be impossible
                attributes[i] = sd.getAttributesMap().get(att);
            }
            attMap.put(att, attributes);
        }

        ScatterPlotData spData = new ScatterPlotData(sampleNames, attMap, dataMap);
        final ScatterPlotFrame igvPlotFrame = new ScatterPlotFrame(" ", spData);
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                igvPlotFrame.setVisible(true);
            }
        });

    }

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
