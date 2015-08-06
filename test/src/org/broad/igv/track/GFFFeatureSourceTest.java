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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012-Jun-22
 */
public class GFFFeatureSourceTest extends AbstractHeadlessTest {

    @Test
    public void testIsGFF() {
        List<String> exts = Arrays.asList("gff3", "gvf", "gff", "gtf");
        for (String ext : exts) {
            String file = "test." + ext;
            assertTrue(GFFFeatureSource.isGFF(file));
            file += ".txt";
            assertTrue(GFFFeatureSource.isGFF(file));
            file += ".gz";
            assertTrue(GFFFeatureSource.isGFF(file));
        }
    }

    private List<Feature> getGeneFeatures(String filepath, String chr, int start, int end) throws Exception {
        TestUtils.createIndex(filepath);

        GFFFeatureSource source = getGffFeatureSource(filepath);

        Iterator<Feature> feats = source.getFeatures(chr, start, end);
        List<Feature> sourceFeats = new ArrayList<Feature>(2);

        while (feats.hasNext()) {
            Feature feat = feats.next();
            assertEquals(chr, feat.getChr());
            sourceFeats.add(feat);
        }
        return sourceFeats;
    }

    private GFFFeatureSource getGffFeatureSource(String filepath) throws IOException, TribbleIndexNotFoundException {
        TribbleFeatureSource fs = TribbleFeatureSource.getFeatureSource(new ResourceLocator(filepath), genome);
        return new GFFFeatureSource(fs);
    }

    @Test
    public void testGetAll() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/gene.sorted.gff3";
        String chr = "chr1";
        int start = 0;
        int end = Integer.MAX_VALUE / 2;
        List<Feature> sourceFeats = getGeneFeatures(filepath, chr, start, end);

        assertEquals(2, sourceFeats.size());


        GFFFeatureSource source = getGffFeatureSource(filepath);

        Iterator<Feature> iter = source.getFeatures(chr, start, end);
        List<Feature> parserFeats = new ArrayList();
        while (iter.hasNext()) parserFeats.add(iter.next());

        assertEquals(parserFeats.size(), sourceFeats.size());

        int sF = 0;
        for (Feature f : parserFeats) {
            BasicFeature sourceFeat = (BasicFeature) sourceFeats.get(sF);
            BasicFeature bf = (BasicFeature) f;
            assertEquals(bf.getExonCount(), sourceFeat.getExonCount());
            assertEquals(bf.getIdentifier(), sourceFeat.getIdentifier());
            sF++;
        }
    }

    @Test
    public void testQuery_01() throws Exception {
        genome = IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome");
        String filepath = org.broad.igv.util.TestUtils.DATA_DIR + "gff/aliased.sorted.gff";

        String chr = "chr5";
        int start = 120960 - 1;
        int end = 125258;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        int geneCount = 0;
        int rnaCount = 0;
        for (Feature feat : features) {


            //Feature feat = features.next();
            assertEquals(chr, feat.getChr());

            BasicFeature bf = (BasicFeature) feat;

            String id = bf.getIdentifier().toLowerCase();
            if (id.contains("gene")) geneCount++;
            if (id.contains("rna")) rnaCount++;

            if ("gene21".equals(id)) {
                assertEquals(0, bf.getExonCount());
                assertEquals("gene", bf.getType());
            }
            if ("rna22".equals(id)) {
                assertEquals(6, bf.getExonCount());
                assertEquals("mRNA", bf.getType());
            }
        }
        assertEquals(2, geneCount);
        assertEquals(2, rnaCount);

    }

    @Test
    public void testPhaseString() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/gene.sorted.gff3";
        String chr = "chr1";
        int start = 0;
        int end = Integer.MAX_VALUE / 2;
        List<Feature> sourceFeats = getGeneFeatures(filepath, chr, start, end);

        String hasPhaseId = "LOC_Os01g01010.1:exon_1";
        int expPhaseNum = 1;
        boolean checkedHasPhase = false;

        BasicFeature bf = (BasicFeature) sourceFeats.get(0);
        for (Exon exon : bf.getExons()) {
            if (exon.getName().equals(hasPhaseId)) {
                checkedHasPhase = true;
                assertEquals(expPhaseNum, exon.getReadingFrame());
            } else {
                assertEquals(-1, exon.getReadingFrame());
            }
        }

        assertTrue(checkedHasPhase);

    }

    /**
     * Test the canonical EDEN sample file as described at http://www.sequenceontology.org/gff3.shtml
     *
     * @throws Exception
     */
    @Test
    public void testEDENSample() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/canonical.eden.sorted.gff3";
        String chr = "chr1";
        int start = 1000 - 1;
        int end = 10000;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        assertEquals(6, features.size());

        /**
         * We split different coding sequences / alternative translations as different features
         */
        int expmRNAFeats = 4;
        int actmRNAFeats = 0;

        List<String> expUniquemRNAIDs = Arrays.asList("mRNA00001", "mRNA00002", "mRNA00003");
        Set<String> mRNAIds = new HashSet<String>(expUniquemRNAIDs.size());

        /*
chr1	.	mRNA	1300	9000	.	+	.	 ID=mRNA00003;Parent=gene00001;Name=EDEN.3
chr1	.	exon	1300	1500	.	+	.	 ID=exon00001;Parent=mRNA00003
chr1	.	exon	3000	3902	.	+	.	 ID=exon00003;Parent=mRNA00001,mRNA00003
chr1	.	CDS	3301	3902	.	+	0	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
chr1	.	CDS	3391	3902	.	+	0	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
chr1	.	CDS	5000	5500	.	+	1	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
chr1	.	CDS	5000	5500	.	+	1	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
chr1	.	CDS	7000	7600	.	+	1	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
chr1	.	CDS	7000	7600	.	+	1	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
         */
        Set<Integer> mRNA3CdStarts = new HashSet<Integer>(Arrays.asList(3300, 3390));

        for (Feature feature : features) {
            BasicFeature bf = (BasicFeature) feature;

            String ident = bf.getIdentifier();
            List<Exon> exons = bf.getExons();
            String type = bf.getType();
            if (type.equals("mRNA")) {
                actmRNAFeats++;
                mRNAIds.add(bf.getIdentifier());
            } else if (type.equals("gene")) {
                assertEquals(bf.getIdentifier(), "gene00001");
                continue;
            } else {
                continue;
            }

            Exon lastExon = exons.get(exons.size() - 1);
            assertEquals(7600, lastExon.getCdEnd());
            assertEquals(lastExon.getCdStart(), lastExon.getStart());
            assertEquals(9000, lastExon.getEnd());

            int midCDSInd = 2;

            if (ident.equals("mRNA00001")) {
                assertEquals(4, bf.getExonCount());
                assertEquals(1201 - 1, exons.get(0).getCdStart());

                Exon secExon = exons.get(1);
                assertWholeExonCoding(secExon);
                assertEquals(3000 - 1, secExon.getStart());
                assertEquals(3902, secExon.getEnd());

            }
            if (ident.equals("mRNA00002")) {
                assertEquals(3, bf.getExonCount());
                assertEquals(1201 - 1, exons.get(0).getCdStart());
                midCDSInd = 1;
            }
            if (ident.equals("mRNA00003")) {
                assertEquals(4, bf.getExonCount());
                boolean passedCdStart = false;
                for (Exon exon : exons) {
                    if (exon.isNonCoding()) {
                        assertEquals("Entire exon is UTR but has coding region: " + exon.getName(), 0, exon.getCodingLength());
                        assertEquals("Entire exon is UTR but has coding region: " + exon.getName(), 0, exon.getCdEnd() - exon.getCdEnd());
                    } else {
                        //There are two coding sequences which differ only in start position
                        if (!passedCdStart) {
                            int cdStart = exon.getCdStart();
                            assertTrue("Exon cdStart not expected, was " + cdStart, mRNA3CdStarts.contains(cdStart));
                            mRNA3CdStarts.remove(cdStart);
                            passedCdStart = true;
                        }
                    }
                    assertTrue(exon.getName().contains("exon0000"));
                }
                assertTrue(passedCdStart);
            }

            Exon midCDS = exons.get(midCDSInd);
            assertEquals(4999, midCDS.getStart());
            assertEquals(5500, midCDS.getEnd());
            assertEquals(midCDS.getStart(), midCDS.getCdStart());
            assertEquals(midCDS.getEnd(), midCDS.getCdEnd());
        }

        assertEquals(0, mRNA3CdStarts.size());

        assertEquals(expmRNAFeats, actmRNAFeats);
        assertEquals(expUniquemRNAIDs.size(), mRNAIds.size());
        for (String expUniquemRNAID : expUniquemRNAIDs) {
            assertTrue("Expected mRNA id not found in file: " + expUniquemRNAID, mRNAIds.contains(expUniquemRNAID));
        }

    }

    /**
     * Test a GFF file which has CDS features, but no features of type "exon"
     *
     * @throws Exception
     */
    @Test
    public void testNoExons() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/NC_009084.gff";
        String chr = "NC_009084.1";
        int start = 0;
        int end = 11302;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        assertEquals(6, features.size());

        for (Feature feat : features) {
            BasicFeature basicFeature = (BasicFeature) feat;
            if (basicFeature.getType().equals("region")) {
                assertEquals("id0", basicFeature.getIdentifier());
            } else if (basicFeature.getType().equals("gene")) {
                assertEquals(1, basicFeature.getExonCount());
                assertTrue(basicFeature.getIdentifier().contains("gene"));
                assertTrue(basicFeature.getName().contains("A1S_"));

                Exon exon = basicFeature.getExons().get(0);
                assertTrue("Exon name incorrect: " + exon.getName(), exon.getName().contains("YP_00"));

                assertWholeExonCoding(exon);
                assertEquals(0, exon.getReadingFrame());
            } else {
                throw new AssertionError("Unknown feature type: " + basicFeature.getType());
            }

        }
    }

    @Test
    public void testtRNA() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/musa_trna.gff3";
        String chr = "chr1";
        int start = 26766;
        int end = 26848;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        assertEquals(2, features.size());

        BasicFeature gene = null, tRNA = null;
        for (Feature feat : features) {
            if (((BasicFeature) feat).getType().equals("gene")) {
                gene = (BasicFeature) feat;
            } else if (((BasicFeature) feat).getType().equals("tRNA")) {
                tRNA = (BasicFeature) feat;
            }
        }

        assertEquals(gene.getIdentifier(), tRNA.getAttributes().get("Parent"));

        assertEquals(1, tRNA.getExonCount());
        Exon exon = tRNA.getExons().get(0);
        assertEquals(tRNA.getIdentifier(), exon.getAttributes().get("Parent"));
        assertWholeExonNonCoding(exon);
    }

    @Test
    public void testMusa_GSMUA_Achr1G00030_001() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/musa_trna.gff3";
        String chr = "chr1";
        int start = 20900;
        int end = 26317;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        assertEquals(3, features.size());

        for (Feature feat : features) {
            BasicFeature bf = (BasicFeature) feat;
            if (bf.getType().equals("mRNA")) {
                assertEquals(8, bf.getExonCount());
                for (Exon exon : bf.getExons()) {
                    assertEquals(bf.getIdentifier(), exon.getAttributes().get("Parent"));
                    assertWholeExonCoding(exon);
                }
            }
        }
    }

    /**
     * Test a simple description of an mRNA containing a 5'utr and cds in GFF-2.  This was submitted with a bug report,
     * the utr was being excluded
     *
     * @throws Exception
     */
    @Test
    public void testUTR() throws Exception {
        String file = TestUtils.DATA_DIR + "gff/utr.gff";

        BufferedReader reader = null;
        Genome genome = null;
        try {
            FeatureSource source =
                    new GFFFeatureSource(TribbleFeatureSource.getFeatureSource(new ResourceLocator(file), genome));

            Iterator<Feature> iter = source.getFeatures("chr1", 0, Integer.MAX_VALUE);
            List<htsjdk.tribble.Feature> features = new ArrayList<Feature>();
            while (iter.hasNext()) {
                features.add(iter.next());
            }

            assertEquals(1, features.size());

            BasicFeature feature = (BasicFeature) features.get(0);
            List<Exon> exons = feature.getExons();
            assertEquals(2, exons.size());

            Exon firstExon = exons.get(0);
            assertTrue(firstExon.isNonCoding());
            assertEquals(11783, firstExon.getStart());
            assertEquals(12157, firstExon.getEnd());
            assertEquals(12157, firstExon.getCdStart());

            Exon secondExon = exons.get(1);
            assertFalse(secondExon.isNonCoding());
            assertEquals(12157, secondExon.getCdStart());
            assertEquals(12157, secondExon.getStart());
            assertEquals(12994, secondExon.getEnd());

        } finally {
            if (reader != null) reader.close();
        }
    }

    /**
     * Test reconstructions of a simple mRNA from a gff3 file.  This test added for a bug caused by using a GFF2
     * codec with a gff3 file.  IGV allows specification of gff3 by either the directive, or file extension.
     *
     * @throws Exception
     */
    @Test
    public void testRNA() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/rna1.gff3";
        BufferedReader reader = null;
        Genome genome = null;
        try {
            GFFFeatureSource source = getGffFeatureSource(filepath);

            Iterator<Feature> iter = source.getFeatures("chr1", 0, Integer.MAX_VALUE);
            int featureCount = 0;
            Feature feature = null;
            while (iter.hasNext()) {
                feature = iter.next();
                featureCount++;
            }
            assertEquals(1, featureCount);

            List<Exon> exons = ((BasicFeature) feature).getExons();
            assertEquals(6, exons.size());


        } finally {
            if (reader != null) reader.close();
        }

    }

    private void assertWholeExonCoding(Exon exon) {
        assertEquals(exon.getCdStart(), exon.getStart());
        assertEquals(exon.getCdEnd(), exon.getEnd());
        assertFalse(exon.isNonCoding());
    }

    private void assertWholeExonNonCoding(Exon exon) {
        assertEquals(exon.getEnd(), exon.getCdStart());
        assertEquals(0, exon.getCodingLength());
        assertTrue(exon.isNonCoding());
    }
}
