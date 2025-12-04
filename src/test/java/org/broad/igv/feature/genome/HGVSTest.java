/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature.genome;

import htsjdk.tribble.NamedFeature;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.tribble.UCSCGeneTableCodec;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class HGVSTest extends AbstractHeadlessTest {

    private static Genome genome;

    @BeforeClass
    public static void setUpGenome() throws Exception {
        String genomePath = TestUtils.DATA_DIR + "genomes/hg38.json";
        genome = GenomeManager.getInstance().loadGenome(genomePath);

        // The tracks aren't actually loaded, so we need to manually load the gene features
        String featureURL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz";
        UCSCGeneTableCodec codec = new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.GENEPRED, genome);
        try (BufferedReader reader = ParsingUtils.openBufferedReader(featureURL)) {
            String line;
            while ((line = reader.readLine()) != null) {
                NamedFeature feature = codec.decode(line);
                FeatureDB.addFeature(feature, genome);
            }
        }
    }


    @Test
    public void testIsValidHGVS() {
        // Note: This validation is permissive - we only check if we can parse positions, not strict HGVS syntax
        assertTrue(HGVS.isValidHGVS("NC_000017.11:g.7579472C>G"));
        assertTrue(HGVS.isValidHGVS("NC_000017.11:g.7579472"));
        assertTrue(HGVS.isValidHGVS("NM_000546.5:c.215C>G"));
        assertTrue(HGVS.isValidHGVS("ENST00000380152.6:c.215C>G"));
        assertTrue(HGVS.isValidHGVS("NR_046018.2:n.100"));
        assertTrue(HGVS.isValidHGVS("NR_046018.2:n.100A>G"));
        assertFalse(HGVS.isValidHGVS("Invalid_HGVS_String"));
        assertFalse(HGVS.isValidHGVS("No_Colon_Or_Type"));
    }

    @Test
    public void testGenomeSearch() throws Exception {

        String hgvs = "NC_000017.11:g.7579472C>G";
        assertTrue(HGVS.isValidHGVS(hgvs));
        SearchCommand.SearchResult result = HGVS.search(hgvs, genome);
        assertEquals("chr17", result.getChr());
        assertEquals(7579471, result.getStart());

        // From UCSC tests
        // chr1	11850845	11867218	NC_000001.11:g.11850846_11867218dup16373	0	+
        hgvs = "NC_000001.11:g.11850846_11867218dup16373";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(11850845, result.getStart());

        // chr1	16782350	17359598	NC_000001.11:g.16782351_17359598del577248	0	+
        hgvs = "NC_000001.11:g.16782351_17359598del577248";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(16782350, result.getStart());

        // chr1	35570363	35656664	NC_000001.11:g.35570364_35656664dup86301	0	+
        hgvs = "NC_000001.11:g.35570364_35656664dup86301";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(35570363, result.getStart());

        hgvs = "NM_000546.6(TP53):c.815T>G";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr17", result.getChr());
        assertEquals(7673804, result.getStart());

    }

    @Test
    public void testSearchInIntrons() throws Exception {

        // Adapted from UCSC tests
        // chr1	11256193	11256194	NM_004958.3(MTOR):c.505-2A>G	0	-
        String hgvs = "ENST00000361445.9:c.505-2A>G";
        assertTrue(HGVS.isValidHGVS(hgvs));
        SearchCommand.SearchResult result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(11256193, result.getStart());

        hgvs = "NM_004958.4(MTOR):c.505-2A>G";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(11256193, result.getStart());


        // chr1	11966984	11966985	NM_000302.3(PLOD1):c.1651-2delA	0	+
        hgvs = "ENST00000196061.5:c.1651-2delA";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(11966984, result.getStart());

        // chr1	40257002	40257003	NM_005857.4(ZMPSTE24):c.-211-1058C>G	0	+
        hgvs = "ENST00000372759.4:c.-211-1058C>G";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(40257002, result.getStart());

//        chr1	40073531	40073535	NM_000310.3(PPT1):c.*526_*529delATCA	0	-
        hgvs = "ENST00000642050.2:c.*526_*529delATCA";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(40073531, result.getStart());
    }

    @Test
    public void testHistoricalRefseq() throws Exception {

        // Adapted from UCSC tests
        // chr1	11256193	11256194	NM_004958.3(MTOR):c.505-2A>G	0	-
        String hgvs = "NM_004958.3(MTOR):c.505-2A>G";
        assertTrue(HGVS.isValidHGVS(hgvs));
        SearchCommand.SearchResult result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(11256193, result.getStart());

        // chr1	11966984	11966985	NM_000302.3(PLOD1):c.1651-2delA	0	+
        hgvs = "NM_000302.3(PLOD1):c.1651-2delA";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(11966984, result.getStart());

    }


    @Test
    public void testProteinHGVS() {
        // Example single codon substitution (hypothetical)
        String hgvs = "ENST00000361445.9:p.77"; // Should map to codon for amino acid 77 of transcript
        assertTrue(HGVS.isValidHGVS(hgvs));
        SearchCommand.SearchResult result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        // Basic sanity: genomic span length should be 3 bases for a single codon
        assertEquals(3, result.getEnd() - result.getStart());

        // Range of codons
        hgvs = "ENST00000361445.9:p.77_78";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        // Span should be >= 6 bases (2 codons) accounting for exon boundaries
        assertTrue(result.getEnd() - result.getStart() >= 6);


        // From UCSC
//    chr1	8361232	8361235	NP_036234.3:p.Thr758Serfs	0	-
        hgvs = "ENST00000400908.7:p.Thr758Serfs";
        assertTrue(HGVS.isValidHGVS(hgvs));
        result = HGVS.search(hgvs, genome);
        assertEquals("chr1", result.getChr());
        assertEquals(3, result.getEnd() - result.getStart());
        assertEquals(8361235, result.getEnd());

        //    chr1	155240657	155240660	NP_001005741.1:p.Leu29Alafs*18	0	-
//    chr1	9262270	9262273	NP_004276.2:p.Val320=	0	+
    }


    @Test
    public void testNonCodingTranscriptNotation() {

//        This HGVS notation represents:
//
//        - **`NR_002196.3`**: RefSeq accession for a non-coding RNA transcript (NR prefix indicates non-coding)
//        - **`(H19)`**: Gene name H19
//                - **`n.`**: Non-coding transcript notation (used for non-coding RNAs like lncRNAs, not mRNA)
//        - **`-7080_-1781`**: Range from position -7080 to -1781 relative to the transcript start
//                - Negative positions indicate nucleotides upstream of the transcript start (position 1)
//        - This refers to a region in the 5' flanking sequence before the transcript begins
//                - **`del`**: Deletion variant
//
//        So this describes a **deletion spanning from 7080 bases upstream to 1781 bases upstream** of the H19
//        non-coding RNA transcript start site. The deletion is 5,300 bases long (7080 - 1781 + 1 = 5,300)
//        in the upstream regulatory region.

        String hgvs = "NR_002196.3(H19):n.-7080_-1781del";
        assertTrue(HGVS.isValidHGVS(hgvs));
        SearchCommand.SearchResult result = HGVS.search(hgvs, genome);
        assertEquals("chr11", result.getChr());
        assertEquals(1999623, result.getStart());  // HGVS is 1-based
        assertEquals(2004923, result.getEnd());


    }

}
