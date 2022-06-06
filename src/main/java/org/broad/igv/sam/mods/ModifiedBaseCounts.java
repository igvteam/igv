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

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 2/6/12
 */
public class ModifiedBaseCounts {


    LinkedHashSet<String> allModifications = new LinkedHashSet<>();
    Map<String, Map<Integer, Integer>> counts;

    Map<String, Map<Integer, Integer>> likelihoods;

    public ModifiedBaseCounts() {
        counts = new HashMap<>();
        likelihoods = new HashMap<>();
    }

    public void incrementCounts(Alignment alignment) {

        // Only works with block formats
        if (alignment.getAlignmentBlocks() == null) return;

        Map<Integer, BaseModification> baseModifications = alignment.getBaseModificationMap();
        if (baseModifications != null) {

            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {
                for (int i = block.getBases().startOffset; i < block.getBases().startOffset + block.getBases().length; i++) {
                    if (baseModifications.containsKey(i)) {
                        BaseModification mod = baseModifications.get(i);
                        double threshold = 255 * PreferencesManager.getPreferences().getAsFloat("SAM.BASEMOD_THRESHOLD");
                        int lh = Byte.toUnsignedInt(mod.likelihood);
                        //if(lh < threshold) continue;

                        int blockIdx = i - block.getBases().startOffset;
                        int position = block.getStart() + blockIdx;   // genomic position
                        Map<Integer, Integer> modCounts = counts.get(mod.modification);
                        if (modCounts == null) {
                            modCounts = new HashMap<>();
                            counts.put(mod.modification, modCounts);
                        }

                        Map<Integer, Integer> modLikelihoods = likelihoods.get(mod.modification);
                        if(modLikelihoods == null) {
                            modLikelihoods = new HashMap<>();
                            likelihoods.put(mod.modification, modLikelihoods);
                        }

                        int c = modCounts.containsKey(position) ? modCounts.get(position) + 1 : 1;
                        int l = modLikelihoods.containsKey(position) ? modLikelihoods.get(position) + lh : lh;

                        modCounts.put(position, c);
                        modLikelihoods.put(position, l);

                        allModifications.add(mod.modification);
                    }
                }
            }
        }
    }

    public int getCount(int position, String modification) {

        Map<Integer, Integer> modCounts = counts.get(modification);
        if (modCounts != null && modCounts.containsKey(position)) {
            return modCounts.get(position);
        } else {
            return 0;
        }
    }

    public int getLikelihood(int position, String modification) {
        Map<Integer, Integer> modLikelihoods = likelihoods.get(modification);
        if (modLikelihoods != null && modLikelihoods.containsKey(position)) {
            return modLikelihoods.get(position);
        } else {
            return getCount(position, modification) * 255;
        }
    }

    public Collection<String> getAllModifications() {
        return allModifications;
    }

    public String getValueString(int position) {
        StringBuffer buffer = new StringBuffer();
        for (Map.Entry<String, Map<Integer, Integer>> entry : counts.entrySet()) {
            String modification = entry.getKey();
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
        for (Map.Entry<String, Map<Integer, Integer>> entry : counts.entrySet()) {

            String modification = entry.getKey();
            Map<Integer, Integer> modCounts = entry.getValue();

            System.out.println("Modification: " + modification);
            for (Map.Entry<Integer, Integer> modKey : modCounts.entrySet()) {
                System.out.println(modKey.getKey() + "  " + modKey.getValue());
            }

        }
    }

}
