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

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.*;
import htsjdk.tribble.index.Index;
import org.broad.igv.Globals;
import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.ChromAliasManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.*;
import org.broad.igv.feature.tribble.reader.IGVComponentMethods;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.IndexCreatorDialog;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.collections.CollUtils;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.vcf.MateVariant;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Jun 27, 2010
 */
abstract public class TribbleFeatureSource implements org.broad.igv.track.FeatureSource {

    IGVFeatureReader reader;
    FeatureReader wrappedReader;
    boolean isVCF;
    Genome genome;

    /**
     * Sequence names used in this file, which might differ from reference sequence names.
     */
    Set<String> sequenceNames;

    ChromAliasManager chromAliasManager;
    private Integer featureWindowSize;
    Object header;
    Class featureClass;

    public static TribbleFeatureSource getFeatureSource(ResourceLocator locator, Genome genome) throws IOException, TribbleIndexNotFoundException {
        return getFeatureSource(locator, genome, true);
    }

    public static TribbleFeatureSource getFeatureSource(ResourceLocator locator, Genome genome, boolean useCache) throws IOException {

        AbstractFeatureReader.setComponentMethods(IGVComponentMethods.getInstance());

        FeatureCodec codec = CodecFactory.getCodec(locator, genome);

        boolean indexExists = false;
        String idxPath = null;
        boolean indexRequired = false;
        // Explicit index path
        if (!locator.isDataURL()) {
            idxPath = locator.getIndexPath();
            if (idxPath != null && idxPath.length() > 0) {
                if (FileUtils.isRemote(idxPath)) {
                    idxPath = HttpUtils.mapURL(idxPath);
                }
                indexExists = true;
            } else {
                String maybeIdxPath = ResourceLocator.indexFile(locator);
                if (maybeIdxPath == null || maybeIdxPath.length() == 0) {
                    indexExists = false;
                } else {
                    if (FileUtils.isRemote(maybeIdxPath)) {
                        maybeIdxPath = HttpUtils.mapURL(maybeIdxPath);
                    }
                    if (FileUtils.resourceExists(maybeIdxPath)) {
                        indexExists = true;
                        idxPath = maybeIdxPath;
                    }
                }
            }

            // Optionally let the user create an index.
            if (!FileUtils.isRemote(locator.getPath())) {
                final int hundredMB = 100000000;
                final int oneGB = 1000000000;
                long size = FileUtils.getLength(locator.getPath());
                indexRequired = size > 10 * oneGB;
                if (!Globals.isHeadless() && locator.isLocal() && !locator.getPath().endsWith(".gz") && !indexExists) {
                    if (size > hundredMB) {
                        createIndex(locator, indexRequired);   // Note, might return null.
                    }
                }
            }
        }

        String path = locator.getPath();
        if (FileUtils.isRemote(path)) {
            path = HttpUtils.createURL(path).toExternalForm();
        }
        if (idxPath != null && idxPath.length() > 0 && FileUtils.isRemote(idxPath)) {
            idxPath = HttpUtils.createURL(idxPath).toExternalForm();
        }
        AbstractFeatureReader basicReader = AbstractFeatureReader.getFeatureReader(path, idxPath, codec, indexRequired || indexExists);


        if (indexExists) { //basicReader.hasIndex()) {
            return new IndexedFeatureSource(basicReader, codec, locator, genome, useCache, true);
        } else {
            return new NonIndexedFeatureSource(basicReader, codec, locator, genome);
        }
    }

    public static AbstractFeatureReader getBasicReader(ResourceLocator locator, Genome genome) throws IOException {
        AbstractFeatureReader.setComponentMethods(IGVComponentMethods.getInstance());
        FeatureCodec codec = CodecFactory.getCodec(locator, genome);

        String path = locator.getPath();
        if (FileUtils.isRemote(path)) {
            path = HttpUtils.mapURL(path);
        }
        return AbstractFeatureReader.getFeatureReader(path, null, codec, false);
    }

    /**
     * Present a dialog for the user to create an index.  This method can return null if the user cancels, or there
     * is an error while creating the index.
     *
     * @param locator
     * @param indexRequired
     * @return
     */
    private static Index createIndex(ResourceLocator locator, boolean indexRequired) {
        File baseFile = new File(locator.getPath());
        File newIdxFile = new File(locator.getPath() + ".idx");
        String messageText = "An index file for " + baseFile.getAbsolutePath() + " could not " +
                "be located. An index is " + (indexRequired ? "required" : "recommended") +
                " to view files of this size.   Click \"Go\" to create one now or \"Cancel" +
                " to proceed without an index.";
        IndexCreatorDialog dialog = IndexCreatorDialog.createShowDialog(IGV.getInstance().getMainFrame(), baseFile, newIdxFile, messageText);
        return (Index) dialog.getIndex();
    }


    private TribbleFeatureSource(AbstractFeatureReader reader, FeatureCodec codec, Genome genome) throws IOException {

        this.genome = genome;
        this.isVCF = codec.getClass() == VCFWrapperCodec.class;
        this.featureClass = codec.getFeatureType();
        this.header = reader.getHeader();
        this.reader = new TribbleReaderWrapper(reader);
        this.wrappedReader = reader;
        this.sequenceNames = new HashSet(reader.getSequenceNames());
        this.chromAliasManager = new ChromAliasManager(sequenceNames, genome);
    }

    protected abstract int estimateFeatureWindowSize(FeatureReader reader);

    protected abstract Collection<String> getSequenceNames();

    public abstract boolean isIndexed();

    public Class getFeatureClass() {
        return featureClass;
    }

    public int getFeatureWindowSize() {
        if (featureWindowSize == null) {
            this.featureWindowSize = estimateFeatureWindowSize(this.wrappedReader);
        }
        return featureWindowSize;
    }

    public Object getHeader() {
        return header;
    }

    static class IndexedFeatureSource extends TribbleFeatureSource {


        private IndexedFeatureSource(AbstractFeatureReader basicReader, FeatureCodec codec, ResourceLocator locator,
                                     Genome genome, boolean useCache, boolean useIndex) throws IOException {
            super(basicReader, codec, genome);
        }


        @Override
        public boolean isIndexed() {
            return true;
        }

        @Override
        public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

            String seqName =  sequenceNames.contains(chr) ? chr : chromAliasManager.getAliasName(chr);
            if (seqName == null) seqName = chr;

            return reader.query(seqName, start, end);
        }

        @Override
        protected Collection<String> getSequenceNames() {
            return reader.getSequenceNames();
        }

        /**
         * Estimate an appropriate feature window size.
         *
         * @param reader
         */
        @Override
        protected int estimateFeatureWindowSize(FeatureReader reader) {

            // Simple formula for VCF.  Appropriate for human 1KG/dbSNp, probably overly conservative otherwise
            if (isVCF) {
                return 10000;
            } else {

            }
            CloseableTribbleIterator<htsjdk.tribble.Feature> iter = null;

            try {
                double mem = RuntimeUtils.getAvailableMemory();
                iter = reader.iterator();
                if (iter.hasNext()) {

                    int nSamples = 1000;
                    htsjdk.tribble.Feature firstFeature = iter.next();
                    htsjdk.tribble.Feature lastFeature = firstFeature;
                    String chr = firstFeature.getChr();
                    int n = 1;
                    long len = 0;
                    while (iter.hasNext() && n < nSamples) {
                        htsjdk.tribble.Feature f = iter.next();
                        if (f != null) {
                            n++;
                            if (f.getChr().equals(chr)) {
                                lastFeature = f;
                            } else {
                                len += lastFeature.getEnd() - firstFeature.getStart() + 1;
                                firstFeature = f;
                                lastFeature = f;
                                chr = f.getChr();
                            }
                        }
                    }
                    double dMem = mem - RuntimeUtils.getAvailableMemory();
                    double bytesPerFeature = Math.max(100, dMem / n);

                    len += lastFeature.getEnd() - firstFeature.getStart() + 1;
                    double featuresPerBase = ((double) n) / len;

                    double targetBinMemory = 20000000;  // 20  mega bytes
                    int maxBinSize = Integer.MAX_VALUE;
                    int bs = Math.min(maxBinSize, (int) (targetBinMemory / (bytesPerFeature * featuresPerBase)));
                    return Math.max(1000000, bs);
                } else {
                    return Integer.MAX_VALUE;
                }
            } catch (IOException e) {
                return 1000000;
            } finally {
                if (iter != null) iter.close();
            }
        }

    }


    public static class NonIndexedFeatureSource extends TribbleFeatureSource {

        /**
         * Map containing all features.  Used only when there is no index.
         */
        Map<String, List<Feature>> featureMap;

        CoverageDataSource coverageData;

        private NonIndexedFeatureSource(AbstractFeatureReader basicReader, FeatureCodec codec, ResourceLocator locator, Genome genome) throws IOException {

            super(basicReader, codec, genome);

            featureMap = new HashMap<>(25);

            try (CloseableIterator<Feature> iter = reader.iterator()) {
                while (iter.hasNext()) {
                    Feature f = iter.next();
                    if (f == null) continue;

                    String seqName = f.getChr();
                    String igvChr = genome == null ? seqName : genome.getCanonicalChrName(seqName);

                    List<Feature> featureList = featureMap.computeIfAbsent(igvChr, k -> new ArrayList<>());
                    featureList.add(f);
                    if (f instanceof IGVNamedFeature named) FeatureDB.addFeature(named, genome);

                    if (this.isVCF && f instanceof Variant v) {
                        String chr2 = v.getAttributeAsString("CHR2");
                        String pos2 = v.getAttributeAsString("END");
                        if (chr2 != null && pos2 != null) {
                            String mateChr = genome == null ? chr2 : genome.getCanonicalChrName(chr2);
                            MateVariant mate = new MateVariant(mateChr, Integer.parseInt(pos2), v);
                            featureList = featureMap.computeIfAbsent(mateChr, k -> new ArrayList<>());
                            featureList.add(mate);
                        }
                    }
                }
            }

            for (List<Feature> featureList : featureMap.values()) {
                FeatureUtils.sortFeatureList(featureList);
            }

            if (genome != null) {
                coverageData = new CoverageDataSource(genome);
                coverageData.computeGenomeCoverage();
                sampleGenomeFeatures();
            }
        }

        @Override
        public boolean isIndexed() {
            return false;
        }

        @Override
        public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
            return coverageData == null ? Collections.<LocusScore>emptyList() :
                    coverageData.getSummaryScoresForRange(chr, start, end, zoom);
        }

        @Override
        public Iterator getFeatures(String chr, int start, int end) throws IOException {

            String seqName =  sequenceNames.contains(chr) ? chr : chromAliasManager.getAliasName(chr);
            if (seqName == null) seqName = chr;


            List<Feature> features = featureMap.get(seqName);
            if (features == null) {
                return Collections.<Feature>emptyList().iterator();
            }
            List<Feature> filteredFeatures = CollUtils.filter(features, FeatureUtils.getOverlapPredicate(seqName, start, end));
            return filteredFeatures.iterator();

        }

        public List<Feature> getAllFeatures() throws IOException {
            List<Feature> allFeatures = new ArrayList<>();
            for (List<Feature> chrFeatures : featureMap.values()) {
                allFeatures.addAll(chrFeatures);
            }
            return allFeatures;
        }


        @Override
        protected Collection<String> getSequenceNames() {
            return featureMap.keySet();
        }

        @Override
        protected int estimateFeatureWindowSize(FeatureReader reader) {
            return -1;
        }

        protected void sampleGenomeFeatures() {
            List<Feature> chrAllFeatures = new ArrayList(1000);
            int sampleLength = (int) ((double) genome.getWGLength() / (1000 * 700));
            int lastFeaturePosition = -1;
            for (String chr : genome.getLongChromosomeNames()) {
                List<Feature> features = featureMap.get(chr);
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

            featureMap.put(Globals.CHR_ALL, chrAllFeatures);
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


                double step = ((double) genome.getWGLength() / 1000) / nBins;
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
}
