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

import org.broad.igv.feature.Range;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;


/**
 * @author jacob
 * @date 2013-Apr-18
 */
public class FpkmTrackingCodecTest {

    @Test
    public void testsample_genesFPKM() throws Exception{
        String path = TestUtils.DATA_DIR + "cufflinks/sample_genes.fpkm_tracking";

        List<? extends Range> values = CufflinksParser.parse(path);

        String[] expGenes = new String[]{"STPG1", "HS3ST1","CFLAR","TFPI", "NDUFAF7"};
        int index = 0;

        for(String expGene: expGenes){
            FPKMValue value = (FPKMValue) values.get(index++);
            assertEquals(expGene, value.getGene());
            assertEquals(2, value.getNumSamples());
        }

        FPKMValue value = (FPKMValue) values.get(0);
        //Sample0
        //0.0585244	0	1.46311	OK
        FPKMSampleValue sampleVal0 = value.getSampleValue(0);
        assertEquals(0.0585244f, sampleVal0.fpkm);

        //Sample1
        //0.01	0.0001	1.0	OK
        FPKMSampleValue sampleVal1 = value.getSampleValue(1);
        assertEquals(0.01f, sampleVal1.fpkm);
        assertEquals(0.0001f, sampleVal1.fpkmLo);
        assertEquals(1.0f, sampleVal1.fpkmHi);

        assertEquals(sampleVal0.getChr(), sampleVal1.getChr());
        assertEquals(sampleVal0.getStart(), sampleVal1.getStart());
        assertEquals(sampleVal0.getEnd(), sampleVal1.getEnd());
    }
}
