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

import org.broad.igv.data.Interval;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Apr-18
 */
public class ExpDiffCodecTest {

    @Test
    public void testsamplegene_expdiff() throws Exception{
        String path = TestUtils.DATA_DIR + "cufflinks/sample_gene_exp.diff";

        List<? extends Interval> values = CufflinksParser.parse(path);

        String[] expGenes = new String[]{"TSPAN6", "TNMD", "DPM1", "SCYL3"};
        int index = 0;

        for(String expGene: expGenes){
            Interval value = values.get(index++);
            assertTrue(value instanceof ExpDiffValue);
            assertEquals(expGene, ((ExpDiffValue) value).getGene());
        }
    }
}
