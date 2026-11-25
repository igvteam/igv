package org.igv.ucsc.twobit;

import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.ByteOrder;

import static org.junit.Assert.*;

public class TwoBitSequenceTest {


    @Test
    public void readSequenceLocal() throws IOException {

        String expectedSequence = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTCTATCTATCTATCTATCTATCTTTTT" +
                "CCCCCCGGGGGGagagagagactc tagcatcctcctacctcacNNacCNctTGGACNCcCaGGGatttcN" +
                "NNcccNNCCNCgN";

        String testFile = TestUtils.DATA_DIR + "twobit/foo.2bit";
        TwoBitSequence reader = new TwoBitSequence(testFile, null, null);

        int start = 5;
        int end = 100;
        String seqString = new String(reader.readSequence("chr1", start, end));
        assertEquals(expectedSequence.substring(start, end), seqString);
    }

    /**
     * The table places an "R" at chr1:11 (1-based) and a "Y" at chr1:50, replacing an "N" and a "T" respectively.
     */
    @Test
    public void readSequenceWithAmbiguityTable() throws IOException {

        String testFile = TestUtils.DATA_DIR + "twobit/foo.2bit";
        String tableFile = TestUtils.DATA_DIR + "twobit/foo.ambiguity.tsv";
        TwoBitSequence reader = new TwoBitSequence(testFile, null, tableFile);

        assertEquals("NNNNNRNNNN", new String(reader.readSequence("chr1", 5, 15)));
        assertEquals("R", new String(reader.readSequence("chr1", 10, 11)));
        assertEquals("NNNNN", new String(reader.readSequence("chr1", 5, 10)));   // end is exclusive
        assertEquals("NNNNN", new String(reader.readSequence("chr1", 11, 16)));
        assertEquals("NNACYCT", new String(reader.readSequence("chr1", 45, 52)));
    }

    @Test
    public void readSequenceRemote() throws IOException {

        String url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit";
        TwoBitSequence reader = new TwoBitSequence(url, null, null);

        // Non-masked no "N" region  chr1:11,830-11,869
        String expectedSeq = "GATTGCCAGCACCGGGTATCATTCACCATTTTTCTTTTCG";
        byte[] seqbytes = reader.readSequence("chr1", 11829, 11869);
        String seq = new String(seqbytes);
        assertEquals(expectedSeq, seq);

        // "N" region  chr1:86-124
        expectedSeq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        seqbytes = reader.readSequence("chr1", 85, 124);
        seq = new String(seqbytes);
        assertEquals(expectedSeq, seq);

        // partially masked region chr1:120,565,295-120,565,335
        expectedSeq = "TATGAACTTTTGTTCGTTGGTgctcagtcctagaccctttt";
        seqbytes = reader.readSequence("chr1", 120565294, 120565335);
        seq = new String(seqbytes);
        assertEquals(expectedSeq, seq);
    }

    /**
     * Seq names :
     * NC_007194.1
     * NC_007195.1
     * NC_007196.1
     * NC_007197.1
     * NC_007198.1
     * NC_007199.1
     * NC_007200.1
     * NC_007201.1
     *
     * @throws IOException
     */
    @Test
    public void twoBitSequenceIndex() throws IOException {

        String url = TestUtils.DATA_DIR + "twobit/GCF_000002655.1.2bit";

        TwoBitIndex index = new TwoBitIndex(url, ByteOrder.LITTLE_ENDIAN, 8);
        long offset = index.searchForOffset("NC_007197.1");

        assertNotNull(offset);
    }

    @Test
    public void twoBitSequenceWithBPTree() throws IOException {

        String url = TestUtils.DATA_DIR + "twobit/GCF_000002655.1.2bit";
        String indexPath = TestUtils.DATA_DIR + "twobit/GCF_000002655.1.2bit.bpt";

        // No index
        TwoBitSequence reader = new TwoBitSequence(url, null, null);
        String expectedSequence = "GCAGGTATCCAAAGCCAGAGGCCTGGTGCTACACGACTGG";
        byte[] seqbytes = reader.readSequence("NC_007194.1", 1644639, 1644679);
        String seq = new String(seqbytes);
        assertEquals(expectedSequence, seq);

        //With index
        reader = new TwoBitSequence(url, indexPath, null);
        seqbytes = reader.readSequence("NC_007194.1", 1644639, 1644679);
        seq = new String(seqbytes);
        assertEquals(expectedSequence, seq);


    }
}