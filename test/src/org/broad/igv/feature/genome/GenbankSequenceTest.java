package org.broad.igv.feature.genome;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;

import static junit.framework.Assert.assertEquals;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/19/12
 * Time: 1:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class GenbankSequenceTest {


    public static final String CHR = "NT_030059";

    @Test
    public void testReadSequence() throws Exception {

        String testFile = TestUtils.DATA_DIR + "gbk/pten_test.gbk";

        //LOCUS       NT_030059             105338 bp    DNA     linear   CON 28-OCT-2010
        BufferedReader reader = new BufferedReader(new FileReader(testFile));
        String nextLine;
        do {
            nextLine = reader.readLine();
        } while (!nextLine.startsWith("ORIGIN"));


        GenbankSequence genbankSequence = new GenbankSequence(CHR, reader);
        reader.close();

        //       61 ttccgaggcg cccgggctcc cggcgcggcg gcggaggggg cgggcaggcc ggcgggcggt
        int start = 60;
        int end = 70;
        String expectedSequence = "ttccgaggcg";
        byte[] seqbytes = genbankSequence.readSequence(CHR, start, end);
        String sequence = new String(seqbytes);
        assertEquals(expectedSequence, sequence);

        // Test end of sequence
        expectedSequence = "tcttgtca";
        end = 105338;
        start = end -  expectedSequence.length();
        seqbytes = genbankSequence.readSequence(CHR, start, end);
        sequence = new String(seqbytes);
        assertEquals(expectedSequence, sequence);


        assertEquals(105338, genbankSequence.getSequenceLenth());
    }
}
