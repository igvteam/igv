package org.broad.igv.feature.genome;

import java.util.Comparator;

/**
 * Comparator for chromosome names.
 */
public class ChromosomeComparator implements Comparator<String> {

    /**
     * @param chr1
     * @param chr2
     * @return
     */
    public int compare(String chr1, String chr2) {

        try {

            // Special rule -- put the mitochondria at the end
            if (chr1.equals("chrM") || chr1.equals("MT")) {
                return 1;
            } else if (chr2.equals("chrM") || chr2.equals("MT")) {
                return -1;
            }

            // Find the first digit
            int idx1 = findDigitIndex(chr1);
            int idx2 = findDigitIndex(chr2);
            if (idx1 == idx2) {
                String alpha1 = idx1 == -1 ? chr1 : chr1.substring(0, idx1);
                String alpha2 = idx2 == -1 ? chr2 : chr2.substring(0, idx2);
                int alphaCmp = alpha1.compareTo(alpha2);
                if (alphaCmp != 0) {
                    return alphaCmp;
                } else {
                    int dig1 = Integer.parseInt(chr1.substring(idx1));
                    int dig2 = Integer.parseInt(chr2.substring(idx2));
                    return dig1 - dig2;
                }
            } else if (idx1 == -1) {
                return 1;

            } else if (idx2 == -1) {
                return -1;
            }
            return idx1 - idx2;
        } catch (Exception numberFormatException) {
            return 0;
        }

    }

    int findDigitIndex(String chr) {

        int n = chr.length() - 1;
        if (!Character.isDigit(chr.charAt(n))) {
            return -1;
        }

        for (int i = n - 1; i > 0; i--) {
            if (!Character.isDigit(chr.charAt(i))) {
                return i + 1;
            }
        }
        return 0;
    }
}
