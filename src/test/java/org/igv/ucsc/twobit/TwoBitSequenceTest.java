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
     * A read starting where a block of "N"s ends, the later blocks of "N"s must still be recognized.
     */
    @Test
    public void readSequenceFromEndOfNBlock() throws IOException {

        String expectedSequence = "ACTCTATCTATCTATCTATCTATCTTTTTCCCCCCGGGGGGagagagagactctagcatcctcctacctcacNNacCNctTGGACNCcCaGGGatttcNNNcccNNCCNCgN";

        TwoBitSequence reader = new TwoBitSequence(TestUtils.DATA_DIR + "twobit/foo.2bit", null, null);
        assertEquals(expectedSequence, new String(reader.readSequence("chr1", 47, 159)));
    }

    /**
     * The table places an "R" at chr1:11 (1-based), replacing an "N", and a "Y" at chr1:50, which is a "T" and is
     * not replaced.
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
        assertEquals("NNACTCT", new String(reader.readSequence("chr1", 45, 52)));
    }

    /**
     * Tables bundled with IGV restore the ambiguity codes of hg38 and hg19.  Other genomes have no table.
     */
    @Test
    public void readSequenceWithBundledAmbiguityTable() throws IOException {

        TwoBitSequence reader = new TwoBitSequence("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit", null, null);
        assertEquals("tgtgNGCGCA", new String(reader.readSequence("chr3", 16902878, 16902888)));
        reader.loadBundledAmbiguityTable("hg38");
        assertEquals("tgtgBGCGCA", new String(reader.readSequence("chr3", 16902878, 16902888)));

        reader = new TwoBitSequence("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit", null, null);
        reader.loadBundledAmbiguityTable("hg19");
        assertEquals("CTCCRRGCTT", new String(reader.readSequence("chr3", 60830758, 60830768)));

        reader = new TwoBitSequence(TestUtils.DATA_DIR + "twobit/foo.2bit", null, null);
        reader.loadBundledAmbiguityTable("foo");
        assertEquals("NNNNNNNNNN", new String(reader.readSequence("chr1", 5, 15)));

        // A table that failed to load falls back to the bundled table, one that loaded is kept
        reader = new TwoBitSequence("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit", null,
                TestUtils.DATA_DIR + "twobit/missing.ambiguity.tsv");
        reader.loadBundledAmbiguityTable("hg38");
        assertEquals("tgtgBGCGCA", new String(reader.readSequence("chr3", 16902878, 16902888)));

        reader = new TwoBitSequence("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit", null,
                TestUtils.DATA_DIR + "twobit/foo.ambiguity.tsv");
        reader.loadBundledAmbiguityTable("hg38");
        assertEquals("tgtgNGCGCA", new String(reader.readSequence("chr3", 16902878, 16902888)));
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