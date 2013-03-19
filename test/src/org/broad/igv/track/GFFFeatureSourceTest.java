/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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

    private List<Feature> getGeneFeatures(String filepath, String chr, int start, int end) throws Exception{
        TestUtils.createIndex(filepath);

        GFFFeatureSource source = new GFFFeatureSource(filepath, genome);

        Iterator<Feature> feats = source.getFeatures(chr, start, end);
        List<Feature> sourceFeats = new ArrayList<Feature>(2);

        while (feats.hasNext()) {
            Feature feat = feats.next();
            assertEquals(chr, feat.getChr());
            sourceFeats.add(feat);
        }
        return sourceFeats;
    }

    @Test
    public void testGetAll() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/gene.sorted.gff3";
        String chr = "chr1";
        int start = 0;
        int end = Integer.MAX_VALUE / 2;
        List<Feature> sourceFeats = getGeneFeatures(filepath, chr, start, end);

        assertEquals(2, sourceFeats.size());


        GFFParser parser = new GFFParser(filepath);
        List<FeatureTrack> tracks = parser.loadTracks(new ResourceLocator(filepath), genome);

        List<Feature> parserFeats = tracks.get(0).getFeatures(chr, start, end);
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
        genome = IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome", true);
        String filepath = org.broad.igv.util.TestUtils.DATA_DIR + "gff/aliased.sorted.gff";

        String chr = "chr5";
        int start = 120960 - 1;
        int end = 125258;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        int geneCount = 0;
        int rnaCount = 0;
        for (Feature feat: features) {


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
    public void testPhaseString() throws Exception{
        String filepath = TestUtils.DATA_DIR + "gff/gene.sorted.gff3";
        String chr = "chr1";
        int start = 0;
        int end = Integer.MAX_VALUE / 2;
        List<Feature> sourceFeats = getGeneFeatures(filepath, chr, start, end);

        String hasPhaseId = "LOC_Os01g01010.1:exon_1";
        int expPhaseNum = 1;
        boolean checkedHasPhase = false;

        BasicFeature bf = (BasicFeature) sourceFeats.get(0);
        for(Exon exon: bf.getExons()){
            if(exon.getName().equals(hasPhaseId)){
                checkedHasPhase = true;
                assertEquals(expPhaseNum, exon.getReadingFrame());
            }else{
                assertEquals(-1, exon.getReadingFrame());
            }
        }

        assertTrue(checkedHasPhase);

    }

    /**
     * Test the canonical EDEN sample file as described at http://www.sequenceontology.org/gff3.shtml
     * @throws Exception
     */
    @Test
    public void testEDENSample() throws Exception{
        String filepath = TestUtils.DATA_DIR + "gff/canonical.eden.sorted.gff3";
        String chr = "ctg123";
        int start = 1000 - 1;
        int end = 10000;

        List<Feature> features = getGeneFeatures(filepath, chr, start, end);
        assertEquals(6, features.size());

        /**
         * We split different coding sequences / alternative translations as different features
         */
        int expmRNAFeats = 4;
        int actmRNAFeats = 0;

        List<String> expUniquemRNAIDs = Arrays.asList("mRNA00001","mRNA00002","mRNA00003");
        Set<String> mRNAIds = new HashSet<String>(expUniquemRNAIDs.size());

        for(Feature feature: features){
            BasicFeature bf = (BasicFeature) feature;

            String ident = bf.getIdentifier();
            List<Exon> exons = bf.getExons();
            if(bf.getType().equals("mRNA")){
                actmRNAFeats++;
                mRNAIds.add(bf.getIdentifier());
            }else{
                continue;
            }

            assertEquals(7600, exons.get(exons.size()-1).getCdEnd());

            if(ident.equals("mRNA00001")){
                assertEquals(4, bf.getExonCount());
                assertEquals(1201 - 1, exons.get(0).getCdStart());
            }if(ident.equals("mRNA00002")){
                assertEquals(3, bf.getExonCount());
                assertEquals(1201 - 1, exons.get(0).getCdStart());
            }if(ident.equals("mRNA00003")){
                assertEquals(4, bf.getExonCount());
            }
        }

        assertEquals(expmRNAFeats, actmRNAFeats);
        assertEquals(expUniquemRNAIDs.size(), mRNAIds.size());
        for(String expUniquemRNAID: expUniquemRNAIDs){
            assertTrue("Expected mRNA id not found in file: " + expUniquemRNAID, mRNAIds.contains(expUniquemRNAID));
        }

    }
}
