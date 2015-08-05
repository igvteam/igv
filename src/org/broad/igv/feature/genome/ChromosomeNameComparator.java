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

package org.broad.igv.feature.genome;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import java.util.Comparator;

/**
 * Comparator for chromosome names. All pure string comparisons are case insensitive.
 * In general, we compare strings lexicographically, but attempt to include numbers
 * if they are in the same position
 * 0. Mitochondria are sorted to the end  (chrM, MT)
 * 1. If BOTH strings contain a number starting at the same location,
 * we first by the leading string, then sort by that number. Examples:
 * a. "chr1" < "chr10" because 1 < 10
 * b. "chrUn_12" > "chr20" because "chrUn" > "chr20"
 * c. "Alpha5" < "gamma1" because Alpha < gamma (numbers in same location are ignored because
 * string comparison takes precedence)
 * 2. Numeric comparisons are performed recursively if numbers found are the same.
 * For example, "scaffold_v2_100" < "scaffold_v2_1000". The first numbers match (2 == 2),
 * but we then compare the trailing strings, and "_100" < "_1000"
 */
public class ChromosomeNameComparator implements Comparator<String> {

    private static ChromosomeNameComparator instance;

    private Table<String, String, Integer> cache = HashBasedTable.create();

    private ChromosomeNameComparator() {
    }

    public static ChromosomeNameComparator get() {
        if (instance == null) {
            instance = new ChromosomeNameComparator();
        }
        return instance;
    }

    private boolean isMito(String chr) {
        return chr.equalsIgnoreCase("chrM") || chr.equalsIgnoreCase("MT");
    }

    public int compare(String chr0, String chr1) {
        if (cache.contains(chr0, chr1)) {
            return cache.get(chr0, chr1);
        }
        int comparison = compareNonCache(chr0, chr1);

        //Just to make sure cache size doesn't go crazy.
        //In general don't expect more than ~50 chromosomes,
        //which would be 50 choose 2 ~= 1250 mappings
        if (cache.size() < 10000) {
            cache.put(chr0, chr1, comparison);
        }
        return comparison;
    }

    public void resetCache() {
        cache.clear();
    }

    public int compareNonCache(String chr0, String chr1) {


        int[] range0 = findDigitRange(chr0);
        int[] range1 = findDigitRange(chr1);

        if (range0 == null || range1 == null || range0[0] != range1[0]) {
            // Special rule -- put the mitochondria at the end
            boolean mito0 = isMito(chr0);
            boolean mito1 = isMito(chr1);
            if (mito0 && !mito1) {
                return +1;
            } else if (!mito0 && mito1) {
                return -1;
            } else if (mito0 && mito1) {
                return 0;
            }

            return chr0.compareToIgnoreCase(chr1);
        } else {
            String alpha1 = chr0.substring(0, range0[0]);
            String alpha2 = chr1.substring(0, range1[0]);
            int alphaCmp = alpha1.compareToIgnoreCase(alpha2);
            if (alphaCmp != 0) {
                return alphaCmp;
            } else {
                long dig1 = 0;
                long dig2 = 0;
                try {
                    dig1 = Long.parseLong(chr0.substring(range0[0], range0[1]));
                    dig2 = Long.parseLong(chr1.substring(range1[0], range1[1]));
                } catch (NumberFormatException e) {
                    // This can occur if numbers are too large for Long.  In this case revert to alpha compare
                    return chr0.compareTo(chr1);
                }
                if (dig1 > dig2) {
                    return 1;
                } else if (dig1 < dig2) {
                    return -1;
                } else {
                    return compare(chr0.substring(range0[1]), chr1.substring(range1[1]));
                }

            }
        }
//        try {
//            // Find the first digit
//            int idx1 = findDigitIndex(chr1);
//            int idx2 = findDigitIndex(chr2);
//            if (idx1 == idx2) {
//                String alpha1 = idx1 == -1 ? chr1 : chr1.substring(0, idx1);
//                String alpha2 = idx2 == -1 ? chr2 : chr2.substring(0, idx2);
//                int alphaCmp = alpha1.compareTo(alpha2);
//                if (alphaCmp != 0) {
//                    return alphaCmp;
//                } else {
//                    int dig1 = Integer.parseInt(chr1.substring(idx1));
//                    int dig2 = Integer.parseInt(chr2.substring(idx2));
//                    return dig1 - dig2;
//                }
//            } else if (idx1 == -1) {
//                return +1;
//            } else if (idx2 == -1) {
//                return -1;
//            }
//            return idx1 - idx2;
//        } catch (Exception numberFormatException) {
//            return 0;
//        }

    }

    /**
     * Return start/end (inclusive/exclusive) locations of first range in string
     * which represent a digit.
     *
     * @param chr
     * @return
     */
    private int[] findDigitRange(String chr) {
        int[] locs = null;
        int loc = 0;
        for (char c : chr.toCharArray()) {
            if (Character.isDigit(c)) {
                if (locs == null) {
                    locs = new int[]{loc, chr.length()};
                }
            } else if (locs != null) {
                locs[1] = loc;
                break;
            }
            loc++;
        }
        return locs;
    }

//    private int findDigitIndex(String chr) {
//
//        int n = chr.length() - 1;
//        if (!Character.isDigit(chr.charAt(n))) {
//            return -1;
//        }
//
//        for (int i = n - 1; i > 0; i--) {
//            if (!Character.isDigit(chr.charAt(i))) {
//                return i + 1;
//            }
//        }
//        return 0;
//    }

}
