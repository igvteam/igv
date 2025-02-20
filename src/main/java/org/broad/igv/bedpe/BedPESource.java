package org.broad.igv.bedpe;

import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.util.Downsampler;
import org.broad.igv.util.FeatureCache;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class BedPESource implements FeatureSource<BedPE> {

    public static final int MAX_WG_COUNT = 10000;
    private FeatureCache<BedPE> featureCache;
    private List<BedPE> allFeatures;
    private List<BedPE> wgFeatures;
    private int visibilityWindow;

    public BedPESource(List<BedPE> allFeatures, Genome genome) {
        this.allFeatures = allFeatures;
        init(allFeatures, genome);
    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {

        if (chr.equals(Globals.CHR_ALL)) {
            return wgFeatures.iterator();
        } else {
            return featureCache.getFeatures(chr,  start,  end).iterator();
        }
    }


    @Override
    public int getFeatureWindowSize() {
        return this.visibilityWindow;
    }

    private void init(List<BedPE> featureList, Genome genome) {

        featureCache = new FeatureCache<>(featureList, 50);

        wgFeatures = createWGFeatures(featureList, genome);

        allFeatures = featureList;

        Map<String, Integer> featureCounts = new HashMap<>();
        for(BedPE f : featureList) {
            int count = featureCounts.containsKey(f.getChr()) ? featureCounts.get(f.getChr()) : 0;
            featureCounts.put(f.getChr(), count + 1);
        }

        // Compute viz window based on feature density
        int vw = Integer.MAX_VALUE;
        for (Map.Entry<String, Integer> entry : featureCounts.entrySet()) {
            String chr = entry.getKey();
            Chromosome chromosome = genome.getChromosome(chr);
            if (chromosome != null) {
                double f = 100000.0 / entry.getValue();
                vw = Math.min(vw, (int) (f * chromosome.getLength()));
            }
        }

        this.visibilityWindow = vw;
    }

    private List<BedPE> createWGFeatures(List<BedPE> features, Genome genome) {

        int size = Math.min(features.size(), MAX_WG_COUNT);
        List<BedPE> wgFeatures = new ArrayList<>(size);

        List<BedPE> sampledFeatures;
        if (features.size() < MAX_WG_COUNT) {
            sampledFeatures = features;
        } else {
            sampledFeatures = downsampleFeatures(features);
        }

        for (BedPE f : sampledFeatures) {

            BedPE wgFeature = new WGFeature(f, genome);

            wgFeatures.add(wgFeature);

        }
        return wgFeatures;
    }

    public static List<BedPE> downsampleFeatures(List<BedPE> features) {

        if (features.isEmpty()) {
            return Collections.EMPTY_LIST;
        }

        BedPE maxScoreFeature = features.stream()
                .max(Comparator.comparing(BedPE::getScore)).get();

        int nBins = maxScoreFeature.getScore() > 0 ? 5 : 1;  // TODO make a function of total # of features & maxCount?
        double binSize = nBins > 1 ? Math.log10(maxScoreFeature.getScore()) / nBins : Integer.MAX_VALUE;

        // Divide features into bins
        List<BedPE>[] binnedFeatures = new List[nBins];
        int counts[] = new int[nBins];
        for (int i = 0; i < nBins; i++) {
            binnedFeatures[i] = new ArrayList<>();
            counts[i] = 0;
        }
        for (BedPE f : features) {
            if (f.isComplement()) continue;
            int bin = f.getScore() <= 0 ? 0 : (int) Math.min(nBins - 1, Math.floor(Math.log10(f.getScore()) / binSize));
            binnedFeatures[bin].add(f);
            counts[bin]++;
        }

        // Add sampled features from each bin
        int featuresPerBin = MAX_WG_COUNT / nBins;
        List<BedPE> sampledFeatures = new ArrayList<>(MAX_WG_COUNT);
        for (int i = 0; i < nBins; i++) {
            List<BedPE> bfs = binnedFeatures[i];
            sampledFeatures.addAll(Arrays.asList(new Downsampler<BedPEFeature>().sample(bfs.toArray(BedPEFeature[]::new), featuresPerBin)));
        }

        // Be sure we keep the maximum feature
        if (maxScoreFeature != null) {
            sampledFeatures.add(maxScoreFeature);
        }

        return sampledFeatures;
    }
}
