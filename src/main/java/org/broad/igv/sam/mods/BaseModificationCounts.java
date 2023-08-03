package org.broad.igv.sam.mods;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.util.collections.ByteArrayList;

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


    Map<BaseModificationKey, Map<Integer, ByteArrayList>> likelihoods;

    transient float lastThreshold;


    public BaseModificationCounts() {
        allModifications = new LinkedHashSet<>();
        likelihoods = new HashMap<>();
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

                    int blockIdx = i - block.getBases().startOffset;
                    int position = block.getStart() + blockIdx;

                    // Loop through base modification sets
                    char canonicalBase = 0;
                    int noModLH = 255;
                    for (BaseModificationSet bmSet : baseModificationSets) {
                        final Map<Integer, Byte> bmSetLikelihoods = bmSet.getLikelihoods();
                        if (bmSetLikelihoods != null && bmSet.containsPosition(i)) {  // TODO or flag == '.' and readbase = bmSet.canonicalbase
                            byte byteLikelihood = bmSetLikelihoods.get(i);
                            BaseModificationKey modKey = BaseModificationKey.getKey(bmSet.getBase(), bmSet.getStrand(), bmSet.getModification());
                            allModifications.add(modKey);
                            pushLikelihood(position, byteLikelihood, modKey);

                            canonicalBase = bmSet.getCanonicalBase();   // Assumed same for all modifications at this position, modificatons on both bases at a position not supported
                            noModLH -= Byte.toUnsignedInt(byteLikelihood);
                        }
                    }

                    if (canonicalBase != 0) {
                        BaseModificationKey noModKey = BaseModificationKey.getKey(canonicalBase, '+', "NONE_" + canonicalBase);
                        allModifications.add(noModKey);
                        pushLikelihood(position, (byte) noModLH, noModKey);
                    }

                }
            }
        }
    }

    private void pushLikelihood(int position, byte byteLikelihood, BaseModificationKey modKey) {
        Map<Integer, ByteArrayList> t = likelihoods.get(modKey);
        if (t == null) {
            t = new HashMap<>();
            likelihoods.put(modKey, t);
        }
        ByteArrayList byteArrayList = t.get(position);
        if (byteArrayList == null) {
            byteArrayList = new ByteArrayList(100);
            t.put(position, byteArrayList);
        }
        byteArrayList.add(byteLikelihood);
    }

    public int getCount(int position, BaseModificationKey key, float threshold) {
        lastThreshold = threshold;
        float scaledThreshold = threshold *  255;
        Map<Integer, ByteArrayList> t = likelihoods.get(key);
        ByteArrayList byteArrayList = t.get(position);
        if (byteArrayList == null) {
            return 0;
        } else {
            int count = 0;
            for (int i = 0; i < byteArrayList.size(); i++) {
                int lh = Byte.toUnsignedInt(byteArrayList.get(i));
                if (lh > scaledThreshold) {
                    count++;
                }
            }
            return count;
        }
    }


    public Set<BaseModificationKey> getAllModificationKeys() {
        return allModifications;
    }

    public String getValueString(int position, AlignmentTrack.ColorOption colorOption) {
        StringBuffer buffer = new StringBuffer();

        //    /**
        //     * Map for capturing modification likelihood "pileup", key is modification identifier,
        //     * value is map of position -> sum of likelihoods for modifications at that position
        //     */
        //    Map<BaseModificationKey, Map<Integer, Integer>> likelihoodSums;


        buffer.append("<br>---------<br>");
        buffer.append("Modifications with likelihood > " + (lastThreshold * 100) + "%");
        for (Map.Entry<BaseModificationKey, Map<Integer, ByteArrayList>> entry : likelihoods.entrySet()) {
            BaseModificationKey key = entry.getKey();
            Map<Integer, ByteArrayList> t = entry.getValue();
            if (t.containsKey(position)) {
                if(key.modification.startsWith("NONE_")) continue;
                int count = this.getCount(position, key, lastThreshold);
                if(count > 0) {
                    String modName = BaseModificationUtils.modificationName(key.modification);
                    buffer.append("<br>&nbsp;&nbsp;" + modName + " (" + key.base + key.strand + "): " + count );
                }
            }
        }

        return buffer.toString();
    }
}
