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

package org.broad.igv.track;

import org.broad.igv.Globals;
import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.collections.CollUtils;
import htsjdk.tribble.Feature;

import java.util.*;

/**
 * Implementation of FeatureSource that wraps a list or map of features for the
 * entire genome.  Instances are typically created by parsing a bed or gff file.  This
 * is a legacy implementation, and does not scale to large feature tracks.
 * <p/>
 * User: jrobinso
 * Date: Jan 31, 2010
 */
public class FeatureCollectionSource implements FeatureSource {

    private TrackType type;

    private Map<String, List<htsjdk.tribble.Feature>> featureMap;

    CoverageDataSource coverageData;

    Genome genome;

    public FeatureCollectionSource(Iterable<? extends Feature> allFeatures, Genome genome) {
        this.genome = genome;
        initFeatures(allFeatures);
        coverageData = new CoverageDataSource(genome);
        coverageData.computeGenomeCoverage();
        sampleGenomeFeatures();
    }

    public Class getFeatureClass() {
        return IGVFeature.class;
    }

    public List<LocusScore> getCoverageScores(String chr, int startLocation, int endLocation, int zoom) {
        return coverageData == null ? Collections.<LocusScore>emptyList() :
                coverageData.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
    }

    /**
     * Return features over the interval as an iterator
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public Iterator<Feature> getFeatures(String chr, int start, int end) {
        return getFeatureList(chr, start, end).iterator();
    }

    public List<Feature> getFeatureList(String chr, int start, int end) {

        List<Feature> features = featureMap.get(chr);
        if (features == null) {
            return Collections.<Feature>emptyList();
        }
        List<Feature> filteredFeatures = CollUtils.filter(features, FeatureUtils.getOverlapPredicate(chr, start, end));
        return filteredFeatures;
    }


    public List<Feature> getFeatures(String chr) {
        return featureMap.get(chr);
    }

    public Set<String> getChrs() {
        return featureMap.keySet();
    }

    public int getFeatureWindowSize() {
        return 0;
    }

    public void setFeatureWindowSize(int size) {
        // ignored
    }


    private void initFeatures(Iterable<? extends Feature> allFeatures) {
        // Separate features by chromosome

            featureMap = new HashMap();
            for (Feature f : allFeatures) {
                List<Feature> fList = featureMap.get(f.getChr());
                if (fList == null) {
                    fList = new ArrayList();
                    featureMap.put(f.getChr(), fList);
                }
                fList.add(f);
            }

            for (List<Feature> featureList : featureMap.values()) {
                FeatureUtils.sortFeatureList(featureList);
            }

            if (featureMap.size() < 100) {
                sampleGenomeFeatures();
            }
    }


    private void setFeatures(String chr, List<Feature> features) {
        FeatureUtils.sortFeatureList(features);
        featureMap.put(chr, features);
    }

    public TrackType getType() {
        return type;
    }

    public void setType(TrackType type) {
        this.type = type;
    }

    protected void sampleGenomeFeatures() {

        if(genome == null) return;

        List<Feature> chrAllFeatures = new ArrayList(1000);
        int sampleLength = (int) ((double) genome.getNominalLength() / (1000 * 700));
        int lastFeaturePosition = -1;
        for (String chr : genome.getLongChromosomeNames()) {
            List<Feature> features = getFeatures(chr);
            if (features != null) {
                long offset = genome.getCumulativeOffset(chr);
                for (Feature feature : features) {
                    if (feature instanceof IGVFeature) {
                        IGVFeature f = (IGVFeature) feature;
                        int genStart = (int) ((offset + f.getStart()) / 1000);
                        int genEnd = (int) ((offset + f.getEnd()) / 1000);
                        if (genEnd > lastFeaturePosition + sampleLength) {
                            BasicFeature f2 = new BasicFeature(Globals.CHR_ALL, genStart, genEnd);
                            if (f instanceof BasicFeature) {
                                BasicFeature bf = (BasicFeature) f;
                                f2.setThickEnd((int) ((offset + bf.getThickEnd()) / 1000));
                                f2.setThickStart((int) ((offset + bf.getThickStart()) / 1000));
                                f2.setName(f.getName());
                            }
                            chrAllFeatures.add(f2);

                            lastFeaturePosition = genEnd;
                        }
                    }
                }
            }
        }

        setFeatures(Globals.CHR_ALL, chrAllFeatures);
    }


    protected void computeGenomeCoverage() {
        int nBins = 1000;
        int[] starts = new int[nBins];
        int[] ends = new int[nBins];
        float[] values = new float[nBins];
        Arrays.fill(values, 0);

        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        double step = ((double) currentGenome.getNominalLength() / 1000) / nBins;
        for (int i = 0; i < nBins; i++) {
            starts[i] = (int) (i * step);
            ends[i] = (int) ((i + 1) * step);
        }


        for (String chr : currentGenome.getLongChromosomeNames()) {
            List<Feature> features = featureMap.get(chr);
            if (features != null) {
                long offset = currentGenome.getCumulativeOffset(chr);
                for (Feature f : features) {
                    int genStart = (int) ((offset + f.getStart()) / 1000);
                    int genEnd = (int) ((offset + f.getEnd()) / 1000);
                    int binStart = (int) (genStart / step);
                    int binEnd = (int) (genEnd / step);
                    for (int i = binStart; i <= binEnd; i++) {
                        values[i] = values[i] + 1;
                    }

                }
            }
        }

        coverageData.dataCache.put(Globals.CHR_ALL, new DataTile(starts, ends, values, null));
    }

    class CoverageDataSource extends AbstractDataSource {

        int windowSize = 1000;
        double dataMin = 0;
        double dataMax = 0;

        Map<String, DataTile> dataCache = new HashMap();

        CoverageDataSource(Genome genome) {
            super(genome);
        }

        protected DataTile getRawData(String chr, int startLocation, int endLocation) {


            DataTile coverageData = dataCache.get(chr);
            if (coverageData == null) {
                coverageData = computeCoverage(chr, startLocation, endLocation);
                dataCache.put(chr, coverageData);
            }
            return coverageData;

        }

        @Override
        protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public int getLongestFeature(String chr) {
            return windowSize;
        }

        public double getDataMax() {
            return dataMax;
        }

        public double getDataMin() {
            return dataMin;
        }

        public TrackType getTrackType() {
            return TrackType.OTHER;  //To change body of implemented methods use File | Settings | File Templates.
        }

        // This won't work for large track!

        private DataTile computeCoverage(String chr, int start, int end) {

            int nBins = (end - start) / windowSize + 1;
            int[] starts = new int[nBins];
            int[] ends = new int[nBins];
            for (int i = 0; i < nBins; i++) {
                starts[i] = start + i * windowSize;
                ends[i] = starts[i] + windowSize;
            }
            float[] values = new float[nBins];
            List<Feature> features = featureMap.get(chr);
            if (features != null) {
                for (Feature f : features) {
                    int startBin = f.getStart() / windowSize;
                    int endBin = f.getEnd() / windowSize;
                    for (int i = startBin; i < endBin; i++) {
                        values[i] = values[i] + 1;
                        dataMax = Math.max(dataMax, values[i]);
                    }
                }
            }
            return new DataTile(starts, ends, values, null);

        }

        protected void computeGenomeCoverage() {
            int nBins = 1000;
            int[] starts = new int[nBins];
            int[] ends = new int[nBins];
            float[] values = new float[nBins];
            Arrays.fill(values, 0);


            double step = ((double) genome.getNominalLength() / 1000) / nBins;
            for (int i = 0; i < nBins; i++) {
                starts[i] = (int) (i * step);
                ends[i] = (int) ((i + 1) * step);
            }


            for (String chr : genome.getLongChromosomeNames()) {
                List<Feature> features = featureMap.get(chr);
                if (features != null) {
                    long offset = genome.getCumulativeOffset(chr);
                    for (Feature f : features) {
                        int genStart = (int) ((offset + f.getStart()) / 1000);
                        int genEnd = (int) ((offset + f.getEnd()) / 1000);
                        int binStart = Math.min(values.length - 1, (int) (genStart / step));
                        int binEnd = Math.min(values.length - 1, (int) (genEnd / step));
                        for (int i = binStart; i <= binEnd; i++) {
                            values[i] = values[i] + 1;
                            dataMax = Math.max(dataMax, values[i]);
                        }

                    }
                }
            }

            dataCache.put(Globals.CHR_ALL, new DataTile(starts, ends, values, null));
        }


        public String getValueString(String chr, double position, ReferenceFrame frame) {

            int zoom = Math.max(0, frame.getZoom());
            List<LocusScore> scores = getSummaryScoresForRange(chr, (int) position - 10, (int) position + 10, zoom);

            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int minWidth = (int) (2 * bpPerPixel);    /* * */

            if (scores == null) {
                return "";
            } else {
                LocusScore score = (LocusScore) FeatureUtils.getFeatureAt(position, minWidth, scores);
                return score == null ? "" : "Mean count: " + score.getScore();
            }
        }


    }


}
