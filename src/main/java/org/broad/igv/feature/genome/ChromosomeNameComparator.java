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
 * Comparator for chromosome names.
 * * Names starting with "chr" are sorted to top (UCSC)
 * * Names that are numeric, after removing "chr", are sorted above names that are not.  Allowance is made for
 * *     suffixes to numeric chromosomes (e.g. chr2a, chr2b).
 * * The "mito" chromosome is sorted to top of non "chr" chromosomes, i.e. the last of the "chr" chromosomes.
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

    public int compare(String chr1, String chr2) {

        boolean t1 = chr1.toLowerCase().startsWith("chr");
        boolean t2 = chr2.toLowerCase().startsWith("chr");

        if (t1 && !t2) {
            return -1;
        } else if (t2 && !t1) {
            return 1;
        }

        String c1 = chr1;
        String c2 = chr2;
        if (t1 && t2) {
            c1 = chr1.substring(3);
            c2 = chr2.substring(3);
        }
        IntegerWithSuffix intSuffix1 = foo(c1);
        IntegerWithSuffix intSuffix2 = foo(c2);

        // Sort integers above others (including integers with suffixes)
        if (intSuffix1 != null && intSuffix2 == null) {
            return -1;
        } else if (intSuffix2 != null && intSuffix1 == null) {
            return 1;
        }

        // If both are integers sort accoringly
        if (intSuffix1 != null && intSuffix2 != null) {
            if (intSuffix1.i == intSuffix2.i) {
                return intSuffix1.suffix.compareTo(intSuffix2.suffix);
            } else {
                return intSuffix1.i - intSuffix2.i;
            }
        }

        // Neither are integers
        else {
            boolean privledged1 = (t1 || intSuffix1 != null || isMito(chr1)) && c1.length() <= 3;
            boolean privledged2 = (t2 || intSuffix2 != null || isMito(chr2)) && c2.length() <= 3;
            if(privledged1 && !privledged2) {
                return -1;
            }
            if(privledged2 && !privledged1) {
                return 1;
            }


            //"1", "2a", "2b", "10", "MT", "scaffold"
            // {"chr1", "chr2a", "chr2b", "chr10", "chr12", "chrLongName", "chrLongName1", "chrX", "chrM", "scaffold"};
            if (isMito(chr1)) {
                return 1;
            } else if (isMito(chr2)) {
                return -1;
            } else {
                return chr1.compareToIgnoreCase(chr2);
            }
        }

    }


    private static boolean isInteger(String str) {
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (c < '0' || c > '9') return false;
        }
        return true;
    }

    IntegerWithSuffix foo(String str) {
        String intSting = "";
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (c < '0' || c > '9') break;
            intSting += c;
        }
        if (intSting.length() == 0) {
            return null;
        } else {
            int i = Integer.parseInt(intSting);
            String suffix = str.substring(intSting.length());
            return new IntegerWithSuffix(i, suffix);
        }
    }

    static class IntegerWithSuffix {

        public IntegerWithSuffix(int i, String suffix) {
            this.i = i;
            this.suffix = suffix;
        }

        int i;
        String suffix;
    }
}
