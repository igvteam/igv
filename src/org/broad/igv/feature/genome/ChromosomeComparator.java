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
 * sometimes many smaller contigs. "important" is determined
 * by size, set in the constructor, or else by the name.
 * <p/>
 * If one chromosome is important and the other is not, it is
 * considered "less" (so it is sorted earlier), otherwise
 * they are sorted by name.
 * User: jacob
 * Date: 2012-Aug-16
 */
public class ChromosomeComparator implements Comparator<Chromosome> {

    /**
     *
     */
    private final int minSizeImportant;

    /**
     * @param minSizeImportant The minimum size to be considered "important" by default.
     *                         Note that a contig might still be considered important if it is named chrXXX
     */
    public ChromosomeComparator(int minSizeImportant) {
        this.minSizeImportant = minSizeImportant;

    }

    @Override
    public int compare(Chromosome o1, Chromosome o2) {
        boolean o1import = isImportant(o1);
        boolean o2import = isImportant(o2);
        boolean checkNames = (o1import == o2import);

        if (checkNames) {
            return ChromosomeNameComparator.get().compare(o1.getName(), o2.getName());
        } else if (o1import) {
            return -1;
        } else {
            return +1;
        }

    }

    private boolean isImportant(Chromosome chromo) {
        if (chromo.getLength() > minSizeImportant) return true;
        if (chromo.getName().toLowerCase().startsWith("chr") &&
                chromo.getName().length() <= 6) return true;
        return false;
    }

    public static LinkedHashMap<String, Chromosome> sortChromosomeList(List<Chromosome> tmpChromos, int minBig,
                                                                       LinkedHashMap<String, Chromosome> chromosomeMap) {
        chromosomeMap.clear();
        Collections.sort(tmpChromos, new ChromosomeComparator(minBig));
        for (int ii = 0; ii < tmpChromos.size(); ii++) {
            Chromosome chromo = tmpChromos.get(ii);
            chromo.setIndex(ii);
            chromosomeMap.put(chromo.getName(), chromo);
        }
        return chromosomeMap;
    }
}
