package org.broad.igv.sam;


import org.broad.igv.feature.genome.Genome;

import java.util.*;


/**
 * Experimental class for phasing alignments for high ploidy regions.
 *
 * @author Jim Robinson
 */

public class HaplotypeUtils {


    private final AlignmentInterval alignmentInterval;
    Genome genome;

    public HaplotypeUtils(AlignmentInterval alignmentInterval, Genome genome) {
        this.alignmentInterval = alignmentInterval;
        this.genome = genome;
    }

    public void clusterAlignments(String chr, int start, int end) {

        AlignmentCounts counts = this.alignmentInterval.getCounts();

        final byte[] reference = genome.getSequence(chr, start, end);

        List<Integer> snpPos = findVariantPositions(start, end, counts, reference);

        KMeans kMeans = new KMeans(reference, start, end, snpPos);

        kMeans.buildMap(this.alignmentInterval.getAlignmentIterator());

        List<KMeans.V> clusters = kMeans.findCentroids(10);

        // Now assign all labels to a centroid

        int n = 0;
        int max = 10;
        while (true) {
            for (String label : kMeans.alignmentMap.keySet()) {

                double min = Double.MAX_VALUE;
                KMeans.V centroid = null;

                for (KMeans.V c : clusters) {
                    double dist = c.distance(label);
                    if (dist < min) {
                        centroid = c;
                        min = dist;
                    }
                }

                if (centroid != null) {
                    centroid.add(label);
                }
            }

            boolean movement = false;
            for (KMeans.V c : clusters) {
                if (c.movement()) {
                    movement = true;
                    break;
                }
            }

            if (movement && n++ < max) {
                for (KMeans.V c : clusters) {
                    c.reset();
                }
            } else {
                break;
            }
        }

        // Now label alignments
        int j = 1;
        for (KMeans.V c : clusters) {

            String label = "" + j++;

            for (String l : c.allLabels) {

                List<Alignment> alignments = kMeans.alignmentMap.get(l);
                for (Alignment a : alignments) {
                    a.setHaplotypeName(label);
                }
            }
        }

    }

    private List<Integer> findVariantPositions(int start, int end, AlignmentCounts counts, byte[] reference) {

        List<Integer> snpPos = new ArrayList<>();

        for (int i = start; i < end; i++) {

            byte ref = reference[i - start];

            float mismatchCount = getMismatchCount(counts, i, ref);

            if (mismatchCount > 0.2f && mismatchCount < 0.8f) {
                snpPos.add(i);
            }
        }
        return snpPos;
    }

    public float getMismatchCount(AlignmentCounts counts, int pos, byte ref) {


        float mismatchQualitySum = 0;


        if (ref < 96) ref += 32;  // a fast "toLowercase"
        for (char c : BaseAlignmentCounts.nucleotides) {
            if (c != ref && c != 'n') {
                mismatchQualitySum += counts.getCount(pos, (byte) c);
            }
        }
        return mismatchQualitySum / counts.getTotalCount(pos);


    }


    static class KMeans {

        private final byte[] reference;
        private final int start;
        private final int end;
        private final List<Integer> positions;
        private final Map<String, List<Alignment>> alignmentMap;

        public KMeans(byte[] reference, int start, int end, List<Integer> positions) {

            this.reference = reference;
            this.start = start;
            this.end = end;
            this.positions = positions;
            this.alignmentMap = new HashMap<>();
        }


        public void buildMap(Iterator<Alignment> iter) {

            while (iter.hasNext()) {

                Alignment alignment = iter.next();

                if (this.start >= alignment.getStart() && this.end <= alignment.getEnd()) {

                    String hapName = "";
                    int dist = 0;

                    for (Integer pos : this.positions) {

                        byte ref = reference[pos - start];
                        boolean found = false;
                        for (AlignmentBlock block : alignment.getAlignmentBlocks()) {

                            if (block.isSoftClipped()) continue;
                            if (block.contains(pos)) {
                                int blockOffset = pos - block.getStart();
                                hapName += (char) block.getBase(blockOffset);
                                found = true;

                                if (ref != block.getBase(blockOffset)) {
                                    dist++;
                                }

                                break;
                            }
                        }
                        if (!found) {
                            hapName += "_";
                        }
                    }

                    hapName = hapName.toLowerCase();

                    List<Alignment> alignments = alignmentMap.get(hapName);
                    if (alignments == null) {
                        alignments = new ArrayList<>();
                        alignmentMap.put(hapName, alignments);
                    }
                    alignments.add(alignment);

                }
            }
        }


        private List<V> findCentroids(int n) {

            List<Map.Entry<String, List<Alignment>>> entries = new ArrayList(alignmentMap.entrySet());

            entries.sort((o1, o2) -> {

                return o2.getValue().size() - o1.getValue().size();

            });


            List<V> classes = new ArrayList<>();
            for (int i = 0; i < n; i++) {

                Map.Entry<String, List<Alignment>> e = entries.get(i);
                V v = new V(e.getKey());
                classes.add(v);
            }

            return classes;
        }


        static class V {

            static byte[] foo = {'a', 'c', 't', 'g', '_'};

            int n;
            int total;
            Map<Byte, int[]> counts;
            byte[] label;

            Set<String> allLabels;
            Set<String> previousLabels;

            public V(String s) {

                this.label = s.toLowerCase().getBytes();

                n = label.length;

                counts = new HashMap<>();
                for (byte b : foo) {
                    counts.put(b, new int[n]);
                }

                allLabels = new HashSet<>();
                previousLabels = new HashSet();

                add(s);
            }

            public void add(String s) {

                byte[] m = s.getBytes();

                if (m.length != n) {
                    System.err.println("Wrong length");
                    return;
                }

                total++;

                for (int i = 0; i < n; i++) {

                    byte b = m[i];
                    if (b < 95) b += 32;  // a fast "toLowercase"

                    int[] cts = counts.get(b);
                    if (cts != null) {
                        cts[i]++;
                    } else {
                        System.err.println("Unknown nuc: " + ((char) m[i]));
                    }
                }

                updateLabel();

                allLabels.add(s);
            }

            void updateLabel() {

                for (int i = 0; i < n; i++) {

                    byte bMax = 0;
                    int cMax = 0;

                    for (byte b : foo) {
                        if (b < 95) b += 32;  // a fast "toLowercase"
                        int[] cts = counts.get(b);
                        if (cts == null) {
                            System.out.println("Null: " + ((char) b));
                        }
                        if (cts[i] > cMax) {
                            cMax = counts.get(b)[i];
                            bMax = b;
                        }
                    }

                    label[i] = bMax;
                }
            }

            double distance(String s) {

                if (s.length() != n) {
                    System.out.println("Unequal lengths");
                    return -Integer.MAX_VALUE;
                }

                // Lowercase
                byte[] b = s.toLowerCase().getBytes();

                double d = 0;
                for (int i = 0; i < b.length; i++) {

                    // Uncomment for Hamming distance
                    //  char c = (char) label[i];
                    //   char c1 = s.charAt(i);
                    //  if(c != c1) {
                    //      d += 1.0;
                    //  }

                    // % of counts matching
                    int[] cts = counts.get(b[i]);
                    double pct = cts[i] / ((double) total);
                    d += 1 - pct;

                }

                return d;
            }

            public void reset() {

                counts = new HashMap<>();
                for (byte b : foo) {
                    counts.put(b, new int[n]);
                }

                previousLabels = new HashSet(allLabels);
                allLabels = new HashSet<>();

                for (int i = 0; i < label.length; i++) {
                    counts.get(label[i])[i] = 1;
                }

                total = 1;

            }

            public boolean movement() {

                if (allLabels.size() != previousLabels.size()) {
                    return true;
                }
                for (String l : previousLabels) {
                    if (!allLabels.contains(l)) return true;
                }
                return false;
            }
        }


    }

}
