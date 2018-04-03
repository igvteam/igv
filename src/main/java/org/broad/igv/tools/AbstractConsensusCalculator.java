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

package org.broad.igv.tools;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.BaseAlignmentCounts;
import org.broad.igv.util.ParsingUtils;

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
            Map<String, String> iupacMap = ParsingUtils.loadIUPACMap();
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
