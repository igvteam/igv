/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.tools;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.BaseAlignmentCounts;
import org.broad.igv.tools.motiffinder.MotifFinderDialog;
import org.broad.igv.tools.motiffinder.MotifFinderSource;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Calculates the consensus sequence based on alignments
 * The rules for defining a consensus sequence are left to
 * descendants.
 *
 * @author jacob
 * @date 2013-Jun-24
 */
public abstract class AbstractConsensusCalculator {

    private static Table<Character, Character, Character> degeneracyTable;

    public char calculateConsensusBase(AlignmentCounts counts, int pos) {
        List<BaseFraction> baseFractions = calculateBaseFractions(counts, pos);
        return calculateConsensusBase(baseFractions);
    }

    /**
     * @param baseFractions Map from nucleotide base -> fraction of total. Should add up to 1f. Allowed
     *                      values are those in {@link org.broad.igv.sam.BaseAlignmentCounts#nucleotides}.
     * @return List of {@code BaseFrequency} objects, in descending order by fraction
     */
    protected abstract char calculateConsensusBase(List<BaseFraction> baseFractions);

    protected final List<BaseFraction> calculateBaseFractions(AlignmentCounts counts, int pos) {

        List<BaseFraction> results = new ArrayList<BaseFraction>(5);
        int totalCount = counts.getTotalCount(pos);
        for (char c : BaseAlignmentCounts.nucleotides) {
            int negCount = counts.getNegCount(pos, (byte) c);
            int posCount = counts.getPosCount(pos, (byte) c);
            int count = negCount + posCount;
            float fraction = ((float) count) / totalCount;
            results.add(new BaseFraction(c, fraction));
        }

        //Sort descending by fraction
        Collections.sort(results, Collections.reverseOrder());
        return results;

    }

    protected final char getDegenerateCode(char base0, char base1) {
        Table<Character, Character, Character> degenTable = getDegeneracyTable();
        return degenTable.get(base0, base1);
    }

    protected static Table<Character, Character, Character> getDegeneracyTable() {
        if (degeneracyTable == null) {
            degeneracyTable = HashBasedTable.create(5, 5);
            //TODO Move the IUPAC loading stuff someplace else, these things shouldn't depend on each other
            Map<String, String> iupacMap = MotifFinderDialog.loadMap(MotifFinderSource.class.getResourceAsStream("resources/iupac_regex_table.txt"));
            for (String s : iupacMap.values()) {
                s = s.replace("[", "").replace("]", "").toLowerCase();
                //System.out.println(s);
                String[] tokens = s.split(",");
                if (tokens.length != 3) continue;
                char a = tokens[1].trim().charAt(0);
                char b = tokens[2].trim().charAt(0);
                char c = tokens[0].trim().charAt(0);
                degeneracyTable.put(a, b, c);
                degeneracyTable.put(b, a, c);
            }
        }
        return degeneracyTable;
    }

    /**
     * Class for storing base/frequency statistics.
     * Note: this class has a natural ordering that is inconsistent with equals.
     */
    public static class BaseFraction implements Comparable {

        public final char base;
        public final float fraction;

        public BaseFraction(char base, float fraction) {
            this.base = base;
            this.fraction = fraction;
        }

        @Override
        public int compareTo(Object o) {
            return (int) Math.signum(this.fraction - ((BaseFraction) o).fraction);
        }
    }
}
