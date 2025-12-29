package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.Range;
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

        List<? extends Range> values = CufflinksParser.parse(path);

        String[] expGenes = new String[]{"TSPAN6", "TNMD", "DPM1", "SCYL3"};
        int index = 0;

        for(String expGene: expGenes){
            Range value = values.get(index++);
            assertTrue(value instanceof ExpDiffValue);
            assertEquals(expGene, ((ExpDiffValue) value).getGene());
        }
    }
}
