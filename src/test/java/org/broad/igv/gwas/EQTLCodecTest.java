package org.broad.igv.gwas;

import org.broad.igv.feature.genome.Genome;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/26/13
 *         Time: 1:36 PM
 */
public class EQTLCodecTest {
    @Test
    public void testDecode() throws Exception {

        Genome genome = null;  // Genome not needed for this test

        EQTLCodec codec = new EQTLCodec(genome);

        String line = "rs4700772\t5\t180341845\tENSG00000168903.7\tBTNL3\t180415845\t129.858921455635\t0\t1.61918925670937e-27";

        EQTLFeature f = codec.decode(line);

        assertEquals("rs4700772", f.getSnp());
        assertEquals("ENSG00000168903.7", f.getGeneId());
        assertEquals("BTNL3", f.getGeneName());
        assertEquals("5", f.getChr());
        assertEquals(180341844, f.getStart());
        assertEquals(180341845, f.getEnd());

    }
}
