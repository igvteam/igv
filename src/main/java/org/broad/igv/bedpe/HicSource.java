package org.broad.igv.bedpe;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.hic.ContactRecord;
import org.broad.igv.hic.HicFile;
import org.broad.igv.hic.NormalizationVector;
import org.broad.igv.hic.Region;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

/**
 * Java port of the JavaScript HicSource.
 * Assumes an instance of HicFile is provided via config as "_hicFile".
 */
public class HicSource implements InteractionSource {

    private final Object genome; // keep generic, expects method getChromosomeName(String)
    private final HicFile hicFile;

    private final int binThreshold = 5;
    private final Map<String, RecordCacheEntry> recordCache = new HashMap<>();

    public HicSource(String path, Genome genome) throws IOException {

        this.hicFile = new HicFile(path, genome);
        this.genome = genome;
    }


    /**
     * Get features for the given region, resolution, and normalization.
     */
    @Override
    public List<BedPE> getFeatures(String chr, int start, int end, double bpPerPixel, String normalization, int maxFeatureCount) throws IOException {

        final int binSize = getBinSize(bpPerPixel);

        List<ContactRecord> records = getRecords(chr, start, end, binSize);

        if (records.isEmpty()) {
            return Collections.emptyList();
        }

        //compute values for thresholding
        List<Float> values = new ArrayList<>();
        int binMin = Integer.MAX_VALUE;
        int binMax = 0;
        for (ContactRecord rec : records) {
            int b1 = rec.bin1();
            int b2 = rec.bin2();
            if (Math.abs(b1 - b2) > this.binThreshold) {
                values.add(rec.counts());
            }
            if (b1 < binMin) binMin = b1;
            if (b2 < binMin) binMin = b2;
            if (b1 > binMax) binMax = b1;
            if (b2 > binMax) binMax = b2;
        }

        if(values.isEmpty()) {
            return Collections.emptyList();
        }

        values.sort(Collections.reverseOrder());
        int t = Math.min(maxFeatureCount, values.size() - 1);
        double threshold = values.get(t);
        double min = values.get((int) Math.floor(t * 0.05));
        double max = values.get((int) Math.floor(t * 0.95));


        List<ContactRecord> thresholdRecords = new ArrayList<>();
        List<ContactRecord> significantRecords = new ArrayList<>();
        for (ContactRecord rec : records) {
            int bin1 = rec.bin1();
            int bin2 = rec.bin2();
            if (Math.abs(bin1 - bin2) <= this.binThreshold) continue;

            float counts = rec.counts();
            if (counts > threshold) {
                significantRecords.add(rec);
            } else if (counts == threshold) {
                thresholdRecords.add(rec);
            }
        }

        // Now add in threshold records
        int nThreshold = maxFeatureCount - significantRecords.size();
        if (nThreshold > 0 && !thresholdRecords.isEmpty()) {
            if (thresholdRecords.size() <= nThreshold) {
                significantRecords.addAll(thresholdRecords);
            } else {
                // Reservoir sampling directly into significantRecords
                // 1. Fill the reservoir with the first nThreshold elements
                significantRecords.addAll(thresholdRecords.subList(0, nThreshold));
                int reservoirStart = significantRecords.size() - nThreshold;

                // 2. Iterate through the rest of the items and replace elements in the reservoir
                for (int i = nThreshold; i < thresholdRecords.size(); i++) {
                    int j = ThreadLocalRandom.current().nextInt(0, i + 1);
                    if (j < nThreshold) {
                        significantRecords.set(reservoirStart + j, thresholdRecords.get(i));
                    }
                }
            }
        }

        double[] normVector = null;
        boolean useNormalization = normalization != null && !"NONE".equals(normalization);
        if (useNormalization) {
            NormalizationVector nv = hicFile.getNormalizationVector(normalization, chr, "BP", binSize);
            if (nv == null) {
                useNormalization = false;
            } else {
                normVector = nv.getValues(binMin, binMax);
                if (normVector == null) {
                    useNormalization = false;
                }
            }
        }

        // Convert contact records to features
        List<BedPE> features = new ArrayList<>();
        String c = getChromosomeNameFromGenome(chr);
        for (ContactRecord rec : significantRecords) {

            int bin1 = rec.bin1();
            int bin2 = rec.bin2();
            float value = rec.counts();
            if (useNormalization) {
                double nvnv = 1;
                try {
                    nvnv = normVector[bin1 - binMin] * normVector[bin2 - binMin];
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
                if (!Double.isNaN(nvnv)) {
                    value /= nvnv;
                } else {
                    continue;
                }
            }

            int start1 = bin1 * binSize;
            int end1 = start1 + binSize;
            int start2 = bin2 * binSize;
            int end2 = start2 + binSize;

            HicFeature f = new HicFeature(c, start1, end1, c, start2, end2, rec.counts(), value);
            int score = (max > min) ? (int) Math.round(200 + Math.min(Math.max((f.getValue() - min) / (max - min), 0), 1) * 600) : 800;
            f.setScore(score);
            features.add(f);
        }

        features.sort(Comparator.comparingInt(a -> a.getStart1()));

        return features;
    }

    private int getBinSize(double bpPerPixel) {
        // choose resolution
        List<Integer> resolutions = hicFile.getBpResolutions();
        int index = 0;
        for (int i = resolutions.size() - 1; i >= 0; i--) {
            if (resolutions.get(i) >= bpPerPixel) {
                index = i;
                break;
            }
        }
        int binSize = resolutions.get(index);
        return binSize;
    }

    public List<String> getNormalizationTypes() {
        return hicFile.getNormalizationTypes();
    }

    @Override
    public boolean hasNormalizationVector(String type, String chr, double bpPerPixel) {
        return hicFile.hasNormalizationVector(type, chr, "BP", getBinSize(bpPerPixel));
    }

    /**
     * Get contact records for the given region and bin size.  We need to query 3 blocks of the matrix
     *   (1) the diagonal block corresponding to the region.  This has captured contacts within the region
     *   (2) the block to the left of the diagonal block.  This captures contacts between the region and the upstream adjacent region
     *   (3) the block to the right of the diagonal block.  This captures contacts between the region and the downstream adjacent region
     * @param chr
     * @param start
     * @param end
     * @param binSize
     * @return
     * @throws IOException
     */
    private List<ContactRecord> getRecords(String chr, int start, int end, int binSize) throws IOException {

        final Region region1 = new Region(chr, start, end);
        List<ContactRecord> records = hicFile.getContactRecords(
                region1,
                region1,
                "BP",
                binSize,
                false
        );

        if(start > 0) {
            Region adjacent = new Region(chr, Math.max(0, start - (end - start)), start);
            List<ContactRecord> adjacentRecords = hicFile.getContactRecords(region1, adjacent, "BP", binSize, false);
            records.addAll(adjacentRecords);
        }

        Region adjacent2 = new Region(chr, end, end + (end - start));
        List<ContactRecord> adjacentRecords2 = hicFile.getContactRecords(region1, adjacent2, "BP", binSize, false);
        records.addAll(adjacentRecords2);

        return records;
    }

    private String getChromosomeNameFromGenome(String chr) {
        try {
            // attempt reflective call genome.getChromosomeName(chr)
            java.lang.reflect.Method m = genome.getClass().getMethod("getChromosomeName", String.class);
            Object res = m.invoke(genome, chr);
            return res != null ? res.toString() : chr;
        } catch (Exception e) {
            return chr;
        }
    }

    public String getNVIString() {
        return hicFile.getNVIString();
    }

    public void setNVIString(String nviString) {
        hicFile.setNVIString(nviString);
    }


    private static class RecordCacheEntry {
        final int start;
        final int end;
        final List<ContactRecord> records;

        RecordCacheEntry(int start, int end, List<ContactRecord> records) {
            this.start = start;
            this.end = end;
            this.records = records;
        }
    }

}