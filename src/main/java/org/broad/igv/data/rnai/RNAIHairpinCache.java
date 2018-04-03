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


package org.broad.igv.data.rnai;

//~--- JDK imports ------------------------------------------------------------

import java.util.*;

/**
 * A singleton class for managing RNAi haripin scores.   The class provides a structure to store
 * hairpin scores keyed by batch ID and gene name, and a method to retrieve them.  Hairpin scores
 * are typically added during parsing of the hairpin file,  and retrieved to render a RNAi track
 * or populate popup text.
 *
 * @author jrobinso
 */
public class RNAIHairpinCache {

    private static RNAIHairpinCache theInstance = new RNAIHairpinCache();

    /**
     * A map for containing hairpin scores keyd by batchId and gene name.
     * The first key (outer map) is batchId.   Inner key is gene name.
     * Each gene can multiple hairpin scores (often 5) which are help in
     * a list.
     */
    Map<String, Map<String, Collection<RNAIHairpinValue>>> hairpinScores;

    private RNAIHairpinCache() {
        hairpinScores = new HashMap();
    }

    /**
     * Method description
     *
     * @return
     */
    public static RNAIHairpinCache getInstance() {
        return theInstance;
    }

    /**
     * Add a haripin score for a particular batch and gene to the hairpin data
     * structure.
     *
     * @param batchId
     * @param geneName
     * @param score
     */
    public void addHairpinScore(String batchId, String geneName, RNAIHairpinValue score) {

        Map<String, Collection<RNAIHairpinValue>> scoresForBatch = hairpinScores.get(batchId);

        if (scoresForBatch == null) {
            scoresForBatch = new HashMap();
            hairpinScores.put(batchId, scoresForBatch);
        }

        Collection<RNAIHairpinValue> scoresForGene = scoresForBatch.get(geneName);

        if (scoresForGene == null) {
            scoresForGene = new TreeSet(new RNAIHairpinNameComparator());
            scoresForBatch.put(geneName, scoresForGene);
        }

        scoresForGene.add(score);
    }

    /**
     * Return the list of scores associated with the specified batch and gene.
     *
     * @param batchId
     * @param geneName
     * @return the scores,  or null if there are no matches.
     */
    public Collection<RNAIHairpinValue> getHairpinScores(String batchId, String geneName) {
        Map<String, Collection<RNAIHairpinValue>> scoresForBatch = hairpinScores.get(batchId);

        if (scoresForBatch == null) {
            return null;
        } else {
            return scoresForBatch.get(geneName);
        }

    }

    private class RNAIHairpinNameComparator implements Comparator<RNAIHairpinValue> {
        public int compare(RNAIHairpinValue hp1, RNAIHairpinValue hp2) {
            if (hp1.getName().length() < hp2.getName().length())
                return -1;

            return hp1.getName().compareTo(hp2.getName());
        }
    }
}
