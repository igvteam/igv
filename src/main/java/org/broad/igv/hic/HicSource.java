package org.broad.igv.hic;

 import org.broad.igv.bedpe.BedPE;
 import org.broad.igv.bedpe.BedPEFeature;
 import org.broad.igv.bedpe.InteractionSource;
 import org.broad.igv.feature.genome.Genome;
 import org.broad.igv.util.collections.LRUCache;


 import java.io.IOException;
 import java.util.*;

 /**
  * Java port of the JavaScript HicSource.
  * Assumes an instance of HicFile is provided via config as "_hicFile".
  */
 public class HicSource implements InteractionSource {

     private final Object genome; // keep generic, expects method getChromosomeName(String)
     private final HicFile hicFile;

     private final int binThreshold;
     private final int percentileThreshold;
     private final Map<String, RecordCacheEntry> recordCache = new HashMap<>();
     private final LRUCache<String, NormalizationVector> normVectorCache = new LRUCache<>(10);

     public HicSource(String path, Map<String, Object> config, Genome genome) throws IOException {

         this.hicFile = HicFile.create(path, config);

         this.genome = genome;
         this.binThreshold = (int) (config.getOrDefault("binThreshold", 5));
         this.percentileThreshold = (int) (config.getOrDefault("percentileThreshold", 80));
     }

     public HicFile getHeader() throws IOException {
         hicFile.init();
         return hicFile;
     }

     /**
      * Get features for the given window.
      */
     public List<BedPE> getFeatures(String chr, int start, int end, double bpPerPixel, String normalization) throws IOException {
         hicFile.init();

         // choose resolution
         List<Integer> resolutions = hicFile.getBpResolutions();
         int index = 0;
         for (int i = resolutions.size() - 1; i >= 0; i--) {
             if (resolutions.get(i) >= bpPerPixel) {
                 index = i;
                 break;
             }
         }
         int selectedIndex = index >= 0 ? index : 0;
         int binSize = resolutions.get(selectedIndex);

         List<ContactRecord> records = getRecords(chr, start, end, binSize);

         int binThresholdLocal = this.binThreshold;

         int nvX1 = 0, nvY1 = 0;
         double[] normVector1 = null, normVector2 = null;
         boolean useNormalization = normalization != null && !"NONE".equals(normalization);
         if (useNormalization) {
             NormalizationVector nv = getNormalizationVector(normalization, chr, binSize);
             if (nv != null) {
                 nvX1 = (int) Math.floor((double) start / binSize);
                 nvY1 = nvX1;
                 normVector1 = nv.getValues(nvX1, (int) Math.ceil((double) end / binSize));
                 normVector2 = nv.getValues(nvY1, (int) Math.ceil((double) end / binSize));
                 if (normVector1 == null || normVector2 == null) {
                     useNormalization = false;
                 }
             } else {
                 useNormalization = false;
             }
         }

         // compute values for thresholding
         List<Double> values = new ArrayList<>();
         for (ContactRecord rec : records) {
             int b1 = rec.bin1();
             int b2 = rec.bin2();
             if (Math.abs(b1 - b2) > binThresholdLocal) {
                 double value = rec.counts();
                 if (useNormalization) {
                     double nvnv = normVector1[b1 - nvX1] * normVector2[b2 - nvY1];
                     if (nvnv != 0 && !Double.isNaN(nvnv)) {
                         value /= nvnv;
                     } else {
                         continue;
                     }
                 }
                 values.add(value);
             }
         }

         double threshold = 0, min = 0, max = 0;
         if (values.size() >= 1000) {
             double[] tmm = percentiles(values, this.percentileThreshold, 50000);
             threshold = tmm[0];
             min = tmm[1];
             max = tmm[2];
         }

         List<BedPE> features = new ArrayList<>();
         for (ContactRecord rec : records) {
             int bin1 = rec.bin1();
             int bin2 = rec.bin2();
             if (Math.abs(bin1 - bin2) <= binThresholdLocal) continue;

             double value = rec.counts();
             if (useNormalization) {
                 double nvnv = normVector1[bin1 - nvX1] * normVector2[bin2 - nvY1];
                 if (nvnv != 0 && !Double.isNaN(nvnv)) {
                     value /= nvnv;
                 } else {
                     continue;
                 }
             }

             if (value < threshold) continue;

             // expect genome has method getChromosomeName(String), reflectively call if necessary
             String c = getChromosomeNameFromGenome(chr);

             int start1 = bin1 * binSize;
             int end1 = start1 + binSize;
             int start2 = bin2 * binSize;
             int end2 = start2 + binSize;

             int fstart = Math.min(start1, start2);
             int fend = Math.max(end1, end2);

             int score = (max > min) ? (int) Math.round(200 + Math.min(Math.max((value - min) / (max - min), 0), 1) * 600) : 800;

             BedPE f = new BedPEFeature(c, start1, end1, c, start2, end2, score);
             features.add(f);
         }

         features.sort(Comparator.comparingInt(a -> a.getStart1()));
         return features;
     }

     public NormalizationVector getNormalizationVector(String normalization, String chr, int binSize) throws IOException {
         String cacheKey = String.format("%s_%s_%d", normalization, chr, binSize);
         if (normVectorCache.containsKey(cacheKey)) {
             return normVectorCache.get(cacheKey);
         }
         NormalizationVector nv = hicFile.getNormalizationVector(normalization, chr, "BP", binSize);
         if (nv != null) normVectorCache.put(cacheKey, nv);
         return nv;
     }

     public List<ContactRecord> getRecords(String chr, int start, int end, int binSize) throws IOException {
         String cacheKey = String.format("%s_%d", chr, binSize);
         RecordCacheEntry cached = recordCache.get(cacheKey);
         if (cached != null) {
             if (start >= cached.start && end <= cached.end) {
                 return cached.records;
             }
         }

         int expandedStart = (int) Math.floor((start + 1) / (double) binSize) * binSize;
         int expandedEnd = (int) Math.ceil((end - 1) / (double) binSize) * binSize;

         List<ContactRecord> records = hicFile.getContactRecords(
                 null,
                 new Region(chr, expandedStart, expandedEnd),
                 new Region(chr, expandedStart, expandedEnd),
                 "BP",
                 binSize,
                 false
         );

         recordCache.clear();
         recordCache.put(cacheKey, new RecordCacheEntry(expandedStart, expandedEnd, records));
         return records;
     }

     public boolean supportsWholeGenome() {
         return false;
     }

     // --- helpers & nested classes ---

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

     private static double[] percentiles(List<Double> array, int p, int maxSize) {
         if (array.isEmpty()) return new double[]{0,0,0};
         List<Double> copy = new ArrayList<>(array);
         Collections.sort(copy);
         int percentileIndex = (int) Math.floor(copy.size() * (p / 100.0));
         int maxSizeIndex = Math.max(0, copy.size() - maxSize);
         int k = Math.max(percentileIndex, maxSizeIndex);
         int i = k + (int) Math.floor((copy.size() - k) * (5.0/100.0));
         int j = k + (int) Math.floor((copy.size() - k) * (95.0/100.0));
         double threshold = copy.get(k);
         double min = copy.get(i);
         double max = copy.get(j);
         return new double[]{threshold, min, max};
     }

     private static class RecordCacheEntry {
         final int start;
         final int end;
         final List<ContactRecord> records;
         RecordCacheEntry(int start, int end, List<ContactRecord> records) {
             this.start = start; this.end = end; this.records = records;
         }
     }

 }