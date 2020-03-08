package org.broad.igv.sam;


import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.util.MessageUtils;

import java.util.*;


/**
 * Experimental class for phasing alignments for high ploidy regions.
 *
 * @author Jim Robinson
 */

public class HaplotypeUtils {

    private static Logger log = Logger.getLogger(HaplotypeUtils.class);

    private final AlignmentInterval alignmentInterval;
    Genome genome;

    public HaplotypeUtils(AlignmentInterval alignmentInterval, Genome genome) {
        this.alignmentInterval = alignmentInterval;
        this.genome = genome;
    }

    public boolean clusterAlignments(String chr, int start, int end, int nClasses) {


        try {
            AlignmentCounts counts = this.alignmentInterval.getCounts();

            final byte[] reference = genome.getSequence(chr, start, end);

            // Find snp positions
            List<Integer> snpPos = findVariantPositions(start, end, counts, reference);

            if (snpPos.size() == 0) {
                MessageUtils.showMessage("No variants in selected range.");
                return false;
            }

            if (snpPos.size() < nClasses - 1) {
                nClasses = snpPos.size() + 1;
                MessageUtils.showMessage("Not enough variants, reducing # of clusters: " + nClasses);
            }


            // Adjust start and end to min and max snp positions, there is no information outside these bounds
            start = snpPos.get(0) - 1;
            end = snpPos.get(snpPos.size() - 1) + 1;

            // Label alignments
            Map<String, List<Alignment>> labelAlignmentMap = labelAlignments(start, end, snpPos, reference, this.alignmentInterval.getAlignmentIterator());
            List<String> labels = new ArrayList(labelAlignmentMap.keySet());
            if (labels.size() < nClasses) {
                MessageUtils.showMessage("Not enough features to create " + nClasses + " classes. Max # of classes = " + labels.size());
                return false;
            }

            // Sort labels (entries) by # of associated alignments
            labels.sort((o1, o2) -> {
                return labelAlignmentMap.get(o2).size() - labelAlignmentMap.get(o1).size();
            });

            // Create initial cluster centroids
            List<V> clusters = new ArrayList<>();
            for (int i = 0; i < nClasses; i++) {
                String label = labels.get(i);
                V v = new V(i + 1, label);
                clusters.add(v);
            }

            // Now assign all labels to a cluster

            int n = 0;
            int max = 50;
            while (true) {
                for (String label : labels) {

                    double min = Double.MAX_VALUE;
                    V centroid = null;

                    for (V c : clusters) {
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
                for (V c : clusters) {
                    if (c.movement()) {
                        movement = true;
                        break;
                    }
                }

                if (movement && n++ < max) {
                    for (V c : clusters) {
                        c.reset();
                    }
                } else {
                    break;
                }
            }

            // Now label alignments by cluster
            for (int i = 0; i < clusters.size(); i++) {
                V c = clusters.get(i);
                String label = "" + c.id;
                for (String l : c.allLabels) {

                    List<Alignment> alignments = labelAlignmentMap.get(l);
                    for (Alignment a : alignments) {
                        a.setHaplotypeName(label);
                    }
                }
            }
            return true;
        } catch (Exception e) {
            log.error("Error clustering alignments", e);
            MessageUtils.showMessage("Error clustering alignments: " + e.getMessage());
            return false;
        }
    }

    private List<Integer> findVariantPositions(int start, int end, AlignmentCounts counts, byte[] reference) {

        List<Integer> snpPos = new ArrayList<>();

        for (int i = start; i < end; i++) {

            byte ref = reference[i - start];

            float mismatchCount = getMismatchCount(counts, i, ref);

            if (mismatchCount > 0.2f) {
                snpPos.add(i);
            }
        }
        return snpPos;
    }

    /**
     * Label the alignments by the base value at each snp position over the region defined by [start, end].
     *
     * @param start
     * @param end
     * @param positions
     * @param reference
     * @param iter
     * @return a map of label -> alignment.
     */
    public Map<String, List<Alignment>> labelAlignments(int start, int end, List<Integer> positions, byte[] reference, Iterator<Alignment> iter) {

        Map<String, List<Alignment>> alignmentMap = new HashMap<>();
        while (iter.hasNext()) {
            Alignment alignment = iter.next();
            if (start >= alignment.getStart() && end <= alignment.getEnd()) {
                String hapName = "";
                for (Integer pos : positions) {
                    boolean found = false;
                    for (AlignmentBlock block : alignment.getAlignmentBlocks()) {
                        if (block.isSoftClipped()) continue;
                        if (block.contains(pos)) {
                            int blockOffset = pos - block.getStart();
                            hapName += (char) block.getBase(blockOffset);
                            found = true;
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
        return alignmentMap;
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

//    List<V> combineClusters(List<V> clusters) {
//
//        List<Pair<Integer, Integer>> combine = new ArrayList<>();
//        for (int i = 0; i < clusters.size(); i++) {
//            for (int j = i + 1; j < clusters.size(); j++) {
//                V c1 = clusters.get(i);
//                V c2 = clusters.get(j);
//                double d = c1.distance(c2);
//                //System.out.println("" + c1.id + " - " + c2.id + "  =  " + d);
//                if (d == 0) {
//                    combine.add(new Pair(i, j));
//                }
//            }
//        }
//    }


    static class V {

        static byte[] foo = {'a', 'c', 't', 'g', '_'};

        int id;
        int n;
        int total;
        Map<Byte, int[]> counts;
        byte[] label;

        Set<String> allLabels;
        Set<String> previousLabels;

        public V(int id, String s) {
            this.id = id;
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

        double distance(V v) {

            byte[] s = v.label;

            if (s.length != n) {
                System.out.println("Unequal lengths");
                return -Integer.MAX_VALUE;
            }

            double d = 0;
            for (int i = 0; i < s.length; i++) {
                byte c = label[i];
                byte c1 = s[i];
                if (c != c1) {
                    d += 1.0;
                }
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
