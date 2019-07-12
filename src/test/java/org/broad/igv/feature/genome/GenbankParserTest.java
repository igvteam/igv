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

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/19/12
 * Time: 1:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class GenbankParserTest {



    @Test
    public void testReadSequence() throws Exception {

        String testFile = TestUtils.DATA_DIR + "gbk/pten_test.gbk";

        GenbankParser genbankParser = new GenbankParser(testFile);
        genbankParser.readFeatures(true);

        String expectedChr = "NT_030059";
        assertEquals(expectedChr, genbankParser.getChr());

        String[] expectedTypes = {"gene", "mRNA", "CDS", "gene", "variation", "variation"};
        int[] expectedStarts = {0, 0, 1032, 82042, -79, 10554};
        int[] expectedEnds = {105338, 105338, 102035, 82643, -78, 10555};

        List<Feature> features = genbankParser.getFeatures();
        for (int i=0; i<expectedTypes.length; i++) {
            BasicFeature bf = (BasicFeature) features.get(i);
            assertEquals(expectedTypes[i], bf.getType());
            assertEquals(expectedStarts[i], bf.getStart());
            assertEquals(expectedEnds[i], bf.getEnd());
        }


        //       61 ttccgaggcg cccgggctcc cggcgcggcg gcggaggggg cgggcaggcc ggcgggcggt
        String chr = genbankParser.getChr();
        int start = 60;
        int end = 70;
        String expectedSequence = "ttccgaggcg";
        byte[] seqbytes = genbankParser.getSequence(chr, start, end);
        String sequence = new String(seqbytes);
        assertEquals(expectedSequence, sequence);

        // Test end of sequence
        expectedSequence = "tcttgtca";
        end = 105338;
        start = end - expectedSequence.length();
        seqbytes = genbankParser.getSequence(chr, start, end);
        sequence = new String(seqbytes);
        assertEquals(expectedSequence, sequence);


        assertEquals(105338, genbankParser.getSequenceLenth());


    }

    /**
     *  variation       -78
     /gene="PTEN"
     /replace="a"
     /replace="g"
     /db_xref="HGMD_null"
     variation       10555
     *
     *      gene            1..105338
     /gene="PTEN"
     /note="Derived by automated computational analysis using
     gene prediction method: BestRefseq."
     /db_xref="GeneID:5728"
     /db_xref="HGNC:9588"
     /db_xref="HPRD:03431"
     /db_xref="MIM:601728"
     mRNA            join(1..1111,30588..30672,62076..62120,67609..67652,
     69576..69814,88681..88822,94416..94582,97457..97681,
     101850..105338)
     /gene="PTEN"
     /product="phosphatase and tensin homolog"
     /exception="unclassified transcription discrepancy"
     /note="Derived by automated computational analysis using
     gene prediction method: BestRefseq."
     /transcript_id="NM_000314.4"
     /db_xref="GI:110224474"
     /db_xref="GeneID:5728"
     /db_xref="HGNC:9588"
     /db_xref="HPRD:03431"
     /db_xref="MIM:601728"
     CDS             join(1033..1111,30588..30672,62076..62120,67609..67652,
     69576..69814,88681..88822,94416..94582,97457..97681,
     101850..102035)
     /gene="PTEN"
     /note="Derived by automated computational analysis using
     gene prediction method: BestRefseq."
     /codon_start=1
     /product="phosphatidylinositol-3,4,5-trisphosphate
     3-phosphatase and dual-specificity protein phosphatase
     PTEN"
     /protein_id="NP_000305.3"
     /db_xref="GI:73765544"
     /db_xref="GeneID:5728"
     /db_xref="HGNC:9588"
     /db_xref="HPRD:03431"
     /db_xref="MIM:601728"
     gene            82043..82643
     /gene="RPL11P3"
     /gene_synonym="RPL11_2_1077"
     /note="Derived by automated computational analysis using
     gene prediction method: Curated Genomic."
     /pseudo
     /db_xref="GeneID:100271271"
     /db_xref="HGNC:36342"
     */
}
