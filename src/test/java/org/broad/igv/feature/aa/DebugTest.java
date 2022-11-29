package org.broad.igv.feature.aa;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.UCSCGeneTableCodec;
import org.broad.igv.util.TestUtils;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class DebugTest {

    public static void main(String [] args) throws IOException {
        testFoo();
    }

    public static void testFoo() throws IOException {

        GenomeManager.getInstance().setCurrentGenome(null);
        Genome genome = TestUtils.loadGenome();

        String expectedSeq = "LRAAK";

        UCSCGeneTableCodec codec = new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.GENEPRED, null);
        String str = "1766\tNM_182679\tchr1\t-\t154830723\t154837903\t154831628\t154835462\t8\t154830723,154832496,154832802,154834450,154834642,154834836,154835383,154837824,\t154832265,154832547,154832894,154834537,154834735,154834910,154835520,154837903,\t0\tGPATCH4\tcmpl\tcmpl\t2,2,0,0,0,1,0,-1,";
        BasicFeature feature = codec.decode(str);

        Exon exon = feature.getExons().get(0);
        Exon prevExon = null;
        Exon nextExon = feature.getExons().get(1);
        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        assertNotNull(aaSeq);

        int aaSeqLen = aaSeq.getSequence().size();
        int start = aaSeqLen - expectedSeq.length();
        for (int i = start; i < aaSeqLen - 1; i++) {
            char exp = expectedSeq.charAt(i - start);
            char actual = aaSeq.getSequence().get(i).getAminoAcid().getSymbol();
            assertEquals("i=" + i, exp, actual);
        }
    }

}
