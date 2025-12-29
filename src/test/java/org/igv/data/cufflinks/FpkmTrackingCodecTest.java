package org.igv.data.cufflinks;

import org.igv.feature.Range;
import org.igv.util.TestUtils;
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
