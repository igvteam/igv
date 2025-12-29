package org.igv.sam.cram;

import htsjdk.samtools.*;
import org.igv.feature.genome.GenomeManager;
import org.junit.Before;
import org.junit.Test;

import java.nio.charset.StandardCharsets;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Created by jrobinso on 5/25/16.
 */

public class IGVReferenceSourceTest {

    public static final String FASTA_URL = "https://igv.org/genomes/data/hg38/hg38.fa";
    public static final String EXPECTED_REFERENCE_BASES = "AAACCCAGGGCAAAGAATCTGGCCCTA"; //bases at 22:27198875-271988902
    public static final String GBK_URL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/NC_012920.1.gbk";
    public static final String GBK_EXPECTED_BASES = "tcatttctctaacagcagtaatattaataattttcatgat";

    @Before
    public void setUp() throws Exception {

    }

    @Test
    public void testGetReferenceBasesByRegion() throws Exception {
        GenomeManager.getInstance().loadGenome(IGVReferenceSourceTest.FASTA_URL);
        IGVReferenceSource refSource = new IGVReferenceSource();
        String expected = EXPECTED_REFERENCE_BASES;
        SAMSequenceRecord rec = new SAMSequenceRecord("chr22", 50818468);
        byte[] bases = refSource.getReferenceBasesByRegion(rec, 27198874, expected.length());
        assertEquals(expected.length(), bases.length);
        assertEquals(expected, new String(bases, StandardCharsets.US_ASCII));

        // Test with a chr alias
        rec = new SAMSequenceRecord("22", 50818468);
        bases = refSource.getReferenceBasesByRegion(rec, 27198874, expected.length());
        assertEquals(expected.length(), bases.length);
        assertEquals(expected, new String(bases, StandardCharsets.US_ASCII));
    }

    @Test
    public void testGenbankReference() throws Exception {
        GenomeManager.getInstance().loadGenome(GBK_URL);
        String expected = GBK_EXPECTED_BASES.toUpperCase(); // Seq at NC_012920:7,279-7,318
        IGVReferenceSource refSource = new IGVReferenceSource();
        SAMSequenceRecord rec = new SAMSequenceRecord("NC_012920", 16569);
        byte[] bases = refSource.getReferenceBasesByRegion(rec, 7278, expected.length());
        assertEquals(expected.length(), bases.length);
        assertEquals(expected, new String(bases, StandardCharsets.US_ASCII));

        // Test with a chr alias
        rec = new SAMSequenceRecord("MT", 16569);
        bases = refSource.getReferenceBasesByRegion(rec, 7278, expected.length());
        assertEquals(expected.length(), bases.length);
        assertEquals(expected, new String(bases, StandardCharsets.US_ASCII));
    }


//    @Test
//    public void testCompressedTiming() throws Exception {
//
//        String fastaURL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa";
//        GenomeManager.getInstance().loadGenome(fastaURL, null);
//        IGVReferenceSource refSource = new IGVReferenceSource();
//        SAMSequenceRecord rec = new SAMSequenceRecord("1", 248956422);
//
//        long t = System.currentTimeMillis();
//        byte[] bases = refSource.getReferenceBases(rec, false);
//        assertEquals(248956422, bases.length);
//        long dt = System.currentTimeMillis() - t;
//
//
//        fastaURL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.gz";
//        GenomeManager.getInstance().loadGenome(fastaURL, null);
//        refSource = new IGVReferenceSource();
//        rec = new SAMSequenceRecord("1", 248956422);
//
//        long t1 = System.currentTimeMillis();
//        bases = refSource.getReferenceBases(rec, false);
//        assertEquals(248956422, bases.length);
//        long dt1 = System.currentTimeMillis() - t1;
//
//        System.out.println(dt + "    " + dt1);
//
//        assertTrue(dt1 < dt);
//
//    }
}
