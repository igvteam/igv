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

package org.broad.igv.tdf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-28
 */
public class TDFDataSourceTest extends AbstractHeadlessTest{

    /**
     * If the old and the new comparator would perform the same,
     * we can still use CHR_ALL.
     * @throws Exception
     */
    @Test
    public void testChrAllSameSortValid() throws Exception{
        String[] chromoArray = new String[]{"chr1", "chr10", "chrM", "chr2", "chr20", "chr3", "chr4"};
        List<String> fileChromos = new ArrayList<String>(Arrays.asList(chromoArray));
        List<String> oldSort = new ArrayList<String>(fileChromos);
        List<String> newSort = new ArrayList<String>(fileChromos);

        TDFDataSource.Pre3Sort(oldSort);
        Collections.sort(newSort, ChromosomeNameComparator.get());

        Assert.assertEquals(oldSort, newSort);

        boolean chrAllGood = TDFDataSource.checkChromoNameOrder(oldSort, genome.getChromosomeNames());
        assertTrue(chrAllGood);
    }
}
