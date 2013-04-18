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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.Mutation;
import org.junit.Test;

import java.util.Iterator;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 10:01 AM
 */
public class MutationDataManagerTest extends AbstractHeadlessTest {

    String path = org.broad.igv.util.TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf";


    @Test
    public void testGetFeatures() throws Exception {

        String sample = "TCGA-02-0055-01A-01W";
        String chr = "chr1";
        int start = 97575730;
        int end = 123482409;

        MutationDataManager mgr = new MutationDataManager(path, genome);
        Iterator<Mutation> mutations =  mgr.getFeatures(sample, chr, start, end);

        int mutationCount = 0;
        while(mutations.hasNext()) {
            Mutation m = mutations.next();
            assertEquals(chr, m.getChr());

            if(m.getStart() >= start && m.getEnd() <= end) {
                mutationCount++;
            }
        }

        // There are 2 mutations for this sample in this interval
        assertEquals(2, mutationCount);

    }


}
