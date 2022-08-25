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

package org.broad.igv.sam.mods;

import org.broad.igv.prefs.PreferencesManager;
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
    LinkedHashSet<Key> allModifications;

    /**
     * Map for counts for each modification (e.g. m, h, etc). Key is base+modification identifier, value is map of position -> count
     */
    Map<Key, Map<Integer, Integer>> counts;

    /**
     * Map for capturing modification likelihood "pileup", key is modification identifier,
     * value is map of position -> sum of likelihoods for modifications at that position
     */
    Map<Key, Map<Integer, Integer>> likelihoodSums;

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
                        Key key = new Key(bmset.getBase(), bmset.getStrand(), bmset.getModification());
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

    public int getCount(int position, Key key) {
        Map<Integer, Integer> modCounts = counts.get(key);
        if (modCounts != null && modCounts.containsKey(position)) {
            return modCounts.get(position);
        } else {
            return 0;
        }
    }

    public int getLikelhoodSum(int position, Key key) {
        Map<Integer, Integer> modLikelihoods = likelihoodSums.get(key);
        if (modLikelihoods != null && modLikelihoods.containsKey(position)) {
            return modLikelihoods.get(position);
        } else {
            return getCount(position, key) * 255;
        }
    }

    public Collection<Key> getAllModifications() {
        return allModifications;
    }

    public String getValueString(int position, AlignmentTrack.ColorOption colorOption) {
        StringBuffer buffer = new StringBuffer();
        for (Map.Entry<Key, Map<Integer, Integer>> entry : counts.entrySet()) {
            String modification = entry.getKey().modification;
            Map<Integer, Integer> modCounts = entry.getValue();
            if (modCounts.containsKey(position)) {
                buffer.append("Modification: " + modification + " (" + modCounts.get(position) + ")<br>");
            }
        }
        return buffer.toString();
    }

    /**
     * For debugging
     */
    public void dump() {
        for (Map.Entry<Key, Map<Integer, Integer>> entry : counts.entrySet()) {

            String modification = entry.getKey().toString();
            Map<Integer, Integer> modCounts = entry.getValue();

            System.out.println("Modification: " + modification);
            for (Map.Entry<Integer, Integer> modKey : modCounts.entrySet()) {
                System.out.println(modKey.getKey() + "  " + modKey.getValue());
            }

        }
    }

    public static class Key {
        char base;
        char strand;
        String modification;

        public Key(char base, char strand, String modification) {
            this.base = base;
            this.strand = strand;
            this.modification = modification;
        }

        public char getBase() {
            return base;
        }

        public char getStrand() {
            return strand;
        }

        public String getModification() {
            return modification;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Key key = (Key) o;
            return base == key.base && strand == key.strand && modification.equals(key.modification);
        }

        @Override
        public int hashCode() {
            return Objects.hash(base, strand, modification);
        }

        @Override
        public String toString() {
            return "" + base + strand + modification;
        }
    }

}
