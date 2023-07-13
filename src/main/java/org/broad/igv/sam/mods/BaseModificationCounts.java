package org.broad.igv.sam.mods;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.AlignmentTrack;

import java.util.*;

/**
 * @author Jim Robinson
 * @date 2/6/12
 */
public class BaseModificationCounts {

    /**
     * Set of all modification seen.
     */
    LinkedHashSet<BaseModificationKey> allModifications;

    /**
     * Map for counts for each modification (e.g. m, h, etc). Key is base+modification identifier, value is map of position -> count
     */
    Map<BaseModificationKey, Map<Integer, Integer>> counts;

    /**
     * Map for capturing modification likelihood "pileup", key is modification identifier,
     * value is map of position -> sum of likelihoods for modifications at that position
     */
    Map<BaseModificationKey, Map<Integer, Integer>> likelihoodSums;

    public BaseModificationCounts() {
        allModifications = new LinkedHashSet<>();
        counts = new HashMap<>();
        likelihoodSums = new HashMap<>();
    }

    /**
     * Increment modification counts for each position spanned by the supplied alignments.  Currently both thresholded
     * and total counts are tallied to support different coloring schemes.
     *
     * @param alignment
     */
    public void incrementCounts(Alignment alignment) {

        // Only works with block formats
        if (alignment.getAlignmentBlocks() == null) return;

        List<BaseModificationSet> baseModificationSets = alignment.getBaseModificationSets();
        if (baseModificationSets != null) {

            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {

                // Loop through read sequence index ("i")
                for (int i = block.getBases().startOffset; i < block.getBases().startOffset + block.getBases().length; i++) {

                    for (BaseModificationSet bmset : baseModificationSets) {

                        //String modification = bmset.getModification();
                        BaseModificationKey key = BaseModificationKey.getKey(bmset.getBase(), bmset.getStrand(), bmset.getModification());
                        Map<Integer, Byte> likelihoods = bmset.getLikelihoods();

                        if (bmset.containsPosition(i)) {

                            int lh = Byte.toUnsignedInt(likelihoods.get(i));

                            Map<Integer, Integer> modCounts = counts.get(key);
                            if (modCounts == null) {
                                modCounts = new HashMap<>();
                                counts.put(key, modCounts);
                            }

                            Map<Integer, Integer> modLikelihoods = likelihoodSums.get(key);
                            if (modLikelihoods == null) {
                                modLikelihoods = new HashMap<>();
                                likelihoodSums.put(key, modLikelihoods);
                            }

                            int blockIdx = i - block.getBases().startOffset;
                            int position = block.getStart() + blockIdx;   // genomic position

                            int c = modCounts.containsKey(position) ? modCounts.get(position) + 1 : 1;


                            int l = modLikelihoods.containsKey(position) ? modLikelihoods.get(position) + lh : lh;
                            modCounts.put(position, c);
                            modLikelihoods.put(position, l);

                        }
                        allModifications.add(key);
                    }
                }
            }
        }
    }

    public int getCount(int position, BaseModificationKey key) {
        Map<Integer, Integer> modCounts = counts.get(key);
        if (modCounts != null && modCounts.containsKey(position)) {
            return modCounts.get(position);
        } else {
            return 0;
        }
    }

    public int getLikelhoodSum(int position, BaseModificationKey key) {
        Map<Integer, Integer> modLikelihoods = likelihoodSums.get(key);
        if (modLikelihoods != null && modLikelihoods.containsKey(position)) {
            return modLikelihoods.get(position);
        } else {
            return getCount(position, key) * 255;
        }
    }

    public Set<BaseModificationKey> getAllModificationKeys() {
        return allModifications;
    }

    public String getValueString(int position, AlignmentTrack.ColorOption colorOption) {
        StringBuffer buffer = new StringBuffer();
        for (Map.Entry<BaseModificationKey, Map<Integer, Integer>> entry : counts.entrySet()) {
            BaseModificationKey key = entry.getKey();
            Map<Integer, Integer> modCounts = entry.getValue();
            if (modCounts.containsKey(position)) {
                final Integer count = modCounts.get(position);
                int lh = (int) (((100.0f / 255) * getLikelhoodSum(position, entry.getKey())) / count);
                String modName = BaseModificationUtils.modificationName(key.modification);
                buffer.append("Modification: " + modName + " (" + key.base + key.strand + ", " + count + "  @ " + lh + "%)<br>");
            }
        }
        return buffer.toString();
    }

    /**
     * For debugging
     */
    public void dump() {
        for (Map.Entry<BaseModificationKey, Map<Integer, Integer>> entry : counts.entrySet()) {

            String modification = entry.getKey().toString();
            Map<Integer, Integer> modCounts = entry.getValue();
            System.out.println("Modification: " + modification);
            for (Map.Entry<Integer, Integer> modKey : modCounts.entrySet()) {
                System.out.println(modKey.getKey() + "  " + modKey.getValue());
            }
        }
    }


}
