package org.igv.sam.mods;

import org.igv.sam.Alignment;
import org.igv.sam.AlignmentBlock;
import org.igv.sam.AlignmentTrack;
import org.igv.util.collections.ByteArrayList;

import java.util.*;
import java.util.stream.Collectors;

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
     * Maxixum likelihood (i.e. maximum of all modifications present) for each position and base moodification key*
     */
    Map<BaseModificationKey, Map<Integer, ByteArrayList>> maxLikelihoods;

    /**
     * Maximum likelihood including no-modification (1 - sum(likelihoods)) for each position and base moodification key*
     */
    Map<BaseModificationKey, Map<Integer, ByteArrayList>> nomodLikelihoods;

    transient float lastThreshold;


    public BaseModificationCounts() {
        allModifications = new LinkedHashSet<>();
        maxLikelihoods = new HashMap<>();
        nomodLikelihoods = new HashMap<>();
    }

    /**
     * Increment modification counts for each position spanned by the supplied alignments.
     *
     * @param alignment
     */
    public void incrementCounts(Alignment alignment) {

        // Only works with block formats
        if (alignment.getAlignmentBlocks() == null) return;

        List<BaseModificationSet> baseModificationSets = alignment.getBaseModificationSets();
        if (baseModificationSets != null) {

            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {

                if(block.isSoftClip()) continue;

                // Loop through read sequence for this block
                for (int blockIdx = 0; blockIdx < block.getBases().length; blockIdx++) {

                    int readIdx = block.getBases().startOffset + blockIdx;
                    int position = block.getStart() + blockIdx;

                    // Loop through base modification sets
                    char canonicalBase = 0;
                    int maxLH = -1;
                    BaseModificationKey maxKey = null;
                    int noModLH = 255;
                    for (BaseModificationSet bmSet : baseModificationSets) {
                        final Map<Integer, Byte> bmSetLikelihoods = bmSet.getLikelihoods();
                        if (bmSetLikelihoods != null && bmSet.containsPosition(readIdx)) {
                            byte byteLikelihood = bmSetLikelihoods.get(readIdx);
                            BaseModificationKey modKey = BaseModificationKey.getKey(bmSet.getBase(), bmSet.getStrand(), bmSet.getModification());
                            allModifications.add(modKey);
                            int lh = Byte.toUnsignedInt(byteLikelihood);
                            noModLH -= lh;
                            if (lh > maxLH) {
                                canonicalBase = bmSet.getCanonicalBase();   // This has to be the same for all modifications at this position
                                maxLH = lh;
                                maxKey = modKey;
                            }
                        }
                    }

                    // Take the modification with highest likelihood, which might be the likelihood of no-modification
                    if (canonicalBase != 0) {
                        BaseModificationKey noModKey = BaseModificationKey.getKey(canonicalBase, '+', "NONE_" + canonicalBase);
                        allModifications.add(noModKey);
                        pushLikelihood(position, (byte) maxLH, maxKey, maxLikelihoods);

                        // 2-color counts, which include no-modification
                        if (noModLH > maxLH) {
                            pushLikelihood(position, (byte) noModLH, noModKey, nomodLikelihoods);
                        } else {
                            pushLikelihood(position, (byte) maxLH, maxKey, nomodLikelihoods);
                        }
                    }
                }
            }
        }
    }

    private void pushLikelihood(int position, byte byteLikelihood, BaseModificationKey modKey, Map<BaseModificationKey, Map<Integer, ByteArrayList>> likelihoods) {
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

    public int getCount(int position, BaseModificationKey key, float threshold, boolean includeNoMods) {

        lastThreshold = threshold;
        float scaledThreshold = threshold * 255;

        Map<Integer, ByteArrayList> t = includeNoMods ? nomodLikelihoods.get(key) : maxLikelihoods.get(key);
        if (t == null) {
            return 0;
        }

        ByteArrayList byteArrayList = t.get(position);
        if (byteArrayList == null) {
            return 0;
        } else {
            int count = 0;
            for (int i = 0; i < byteArrayList.size(); i++) {
                int lh = Byte.toUnsignedInt(byteArrayList.get(i));
                if (lh >= scaledThreshold) {
                    count++;
                }
            }
            return count;
        }
    }

    public int getLikelihoodSum(int position, BaseModificationKey key, float threshold, boolean includeNoMods) {
        lastThreshold = threshold;
        float scaledThreshold = threshold * 255;
        Map<Integer, ByteArrayList> t =includeNoMods ? nomodLikelihoods.get(key) :  maxLikelihoods.get(key);
        ByteArrayList byteArrayList = t.get(position);
        if (byteArrayList == null) {
            return 0;
        } else {
            int count = 0;
            for (int i = 0; i < byteArrayList.size(); i++) {
                int lh = Byte.toUnsignedInt(byteArrayList.get(i));
                if (lh >= scaledThreshold) {
                    count += lh;
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
        StringBuffer nomodBuffer = new StringBuffer();

        //    /**
        //     * Map for capturing modification likelihood "pileup", key is modification identifier,
        //     * value is map of position -> sum of likelihoods for modifications at that position
        //     */
        //    Map<BaseModificationKey, Map<Integer, Integer>> likelihoodSums;


        buffer.append("<br>---------<br>");
        buffer.append("Modifications with likelihood > " + (lastThreshold * 100) + "%");

        final boolean includeNomods = colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR;

        Map<BaseModificationKey, Map<Integer, ByteArrayList>> l = includeNomods ? nomodLikelihoods : maxLikelihoods;
        for (BaseModificationKey key : l.keySet()) {
            Map<Integer, ByteArrayList> t = l.get(key);
            if (t.containsKey(position)) {
                int count = this.getCount(position, key, lastThreshold, includeNomods);
                if (count > 0) {
                    int likelihoodSum = getLikelihoodSum(position, key, lastThreshold, includeNomods);
                    int averageLikelihood = (int) ((((double) likelihoodSum) / count) * .3921568);        //.39 => 100/255
                    String modName = BaseModificationUtils.modificationName(key.modification);
                    if(key.modification.startsWith("NONE_")) {
                        nomodBuffer.append("<br>&nbsp;&nbsp;" + modName + ": " + count + "  @ average likelihood " + averageLikelihood + "%");
                    } else {
                        buffer.append("<br>&nbsp;&nbsp;" + modName + " (" + key.base + key.strand + "): " + count + "  @ average likelihood " + averageLikelihood + "%");
                    }
                }
            }
        }
        buffer.append(nomodBuffer.toString());
        return buffer.toString();
    }
}
