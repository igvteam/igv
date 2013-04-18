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

package org.broad.igv.data.cufflinks;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Apr-18
 */
public class FpkmTrackingCodecTest {

    @Test
    public void testsample_genesFPKM() throws Exception{
        String path = TestUtils.DATA_DIR + "cufflinks/sample_genes.fpkm_tracking";

        List<? extends CufflinksValue> values = CufflinksParser.parse(path);

        String[] expGenes = new String[]{"STPG1", "HS3ST1","CFLAR","TFPI", "NDUFAF7"};
        int index = 0;

        for(String expGene: expGenes){
            CufflinksValue value = values.get(index++);
            assertTrue(value instanceof FPKMValue);
            assertEquals(expGene, value.getGene());
        }
    }
}
