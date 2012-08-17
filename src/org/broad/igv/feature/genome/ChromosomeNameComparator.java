/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.genome;

import java.util.Comparator;

/**
 * Comparator for chromosome names. All pure string comparisons are case insensitive,
 * with some exceptions:
 * 0. Mitochondria are sorted to the end
 * 1. If BOTH strings contain a number starting at the same location,
 * we first by the leading string, then sort by that number. Examples:
 * chr1 < chr10 because 1 < 10
 * chrUn_12 > chr20 because chrUn > chr20
 * Alpha5 < gamma1 because Alpha < gamma (numbers in same location are ignored)
 */
public class ChromosomeNameComparator implements Comparator<String> {

    private static ChromosomeNameComparator instance;

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


    public int compare(String chr1, String chr2) {


        // Special rule -- put the mitochondria at the end
        // Assuming we won't have 2
        if (isMito(chr1)) {
            return +1;
        } else if (isMito(chr2)) {
            return -1;
        }

        int[] range1 = findDigitRange(chr1);
        int[] range2 = findDigitRange(chr2);


        if (range1 == null || range2 == null || range1[0] != range2[0]) {
            return chr1.compareToIgnoreCase(chr2);
        } else {
            String alpha1 = chr1.substring(0, range1[0]);
            String alpha2 = chr2.substring(0, range2[0]);
            int alphaCmp = alpha1.compareToIgnoreCase(alpha2);
            if (alphaCmp != 0) {
                return alphaCmp;
            } else {
                int dig1 = Integer.parseInt(chr1.substring(range1[0], range1[1]));
                int dig2 = Integer.parseInt(chr2.substring(range2[0], range2[1]));
                return dig1 - dig2;
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
