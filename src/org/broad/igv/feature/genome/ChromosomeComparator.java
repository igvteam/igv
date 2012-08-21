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
