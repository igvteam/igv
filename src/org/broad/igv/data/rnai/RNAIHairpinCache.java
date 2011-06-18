/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
