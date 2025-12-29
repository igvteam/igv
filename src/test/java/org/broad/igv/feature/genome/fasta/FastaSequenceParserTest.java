package org.broad.igv.feature.genome.fasta;

import org.broad.igv.feature.genome.InMemorySequence;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 3/26/12
 */
public class FastaSequenceParserTest {

    static InMemorySequence fastaSequence;

    @BeforeClass
    public static void setup() throws IOException {
        String path = "https://igv.genepattern.org/test/fasta/chr22.fa";
        Map<String, byte[]> sequenceMap = FastaSequenceParser.parseFasta(path);
        fastaSequence = new InMemorySequence(sequenceMap);
    }

    @Test
    public void testReadSequence() throws Exception {


        String chr = "chr22";
        int start = 31084127;
        int end = 31084167;

        String expectedSequence = "gccaccatgcctggctagttttttgtatttttagtagaga";


        byte[] seq = fastaSequence.getSequence(chr, start, end);
        String seqString = new String(seq);


        assertEquals(expectedSequence, seqString);
    }

    @Test
    public void testReadEnd() throws Exception {

        String chr = "chr22";
        int chrLen = 51304566;
        int start = chrLen - 10;
        int end = chrLen + 10;
        byte[] bytes = fastaSequence.getSequence(chr, start, end);
        assertEquals(10, bytes.length);

        byte[] expectedSequence = "NNNNNNNNNN".getBytes();

        for (int i = 0; i < 10; i++) {
            assertEquals(expectedSequence[i], bytes[i]);
        }
    }


}


