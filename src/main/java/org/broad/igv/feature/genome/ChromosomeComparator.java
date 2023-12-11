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

import org.broad.igv.feature.Chromosome;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * For comparing chromosomes. We order by "important" vs "unimportant".
 * The idea being there are main chromosomes (chr1...chrX, chrM) and
 * sometimes many smaller contigs.
 * <p>
 * Date: 2012-Aug-16
 */
public class ChromosomeComparator implements Comparator<Chromosome> {


    public ChromosomeComparator() {
    }

    @Override
    public int compare(Chromosome o1, Chromosome o2) {

        // Rules based on UCSC naming convention
        boolean o1import = isMajor(o1.getName());
        boolean o2import = isMajor(o2.getName());

        if (o1import && o2import) {
            return ChromosomeNameComparator.get().compare(o1.getName(), o2.getName());
        } else if (o1import) {
            return -1;
        } else if (o2import) {
            return 1;
        } else if (isMyto(o1.getName())) {
            return -1;
        } else if (isMyto(o2.getName())) {
            return 1;
        } else if (o1.getName().startsWith("chr") && !o2.getName().startsWith("chr")) {
            return -1;
        } else if (o2.getName().startsWith("chr") && !o1.getName().startsWith("chr")) {
            return 1;
        } else {
            return o2.getLength() - o1.getLength();
        }
    }

    // Rule based on common UCSC & NCBI naming conventions.  Meant to separate major chromosomes from contings
    boolean isMajor(String chr) {
        if (chr.startsWith("chr")) chr = chr.substring(3);
        return (chr.length() < 3 && isInteger(chr)) || chr.equals("X") || chr.equals("Y");
    }

    boolean isMyto(String chr) {
        return chr.equals("chrM") || chr.equals("MT");
    }

    boolean isInteger(String chr) {
        for (int i = 0; i < chr.length(); i++) {
            if (chr.charAt(i) < '0' || chr.charAt(i) > '9') return false;
        }
        return true;
    }

}
