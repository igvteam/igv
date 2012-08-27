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
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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

    @Test
    public void testGetAll() throws Exception {
        String filepath = org.broad.igv.util.TestUtils.DATA_DIR + "gff/gene.sorted.gff3";
        TestUtils.createIndex(filepath);

        GFFFeatureSource source = new GFFFeatureSource(filepath, genome);
        String chr = "chr1";
        int start = 0;
        int end = Integer.MAX_VALUE / 2;

        Iterator<Feature> feats = source.getFeatures(chr, start, end);
        List<Feature> sourceFeats = new ArrayList<Feature>(2);

        while (feats.hasNext()) {
            Feature feat = feats.next();
            assertEquals(chr, feat.getChr());
            sourceFeats.add(feat);
        }
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
        String filepath = org.broad.igv.util.TestUtils.DATA_DIR + "gff/aliased.sorted.gff";
        TestUtils.createIndex(filepath);
        Genome genome = IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome", true);

        GFFFeatureSource source = new GFFFeatureSource(filepath, genome);
        String chr = "chr5";
        int start = 120960 - 1;
        int end = 125258;

        CloseableTribbleIterator<Feature> features = source.getFeatures(chr, start, end);
        int geneCount = 0;
        int rnaCount = 0;
        while (features.hasNext()) {


            Feature feat = features.next();
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
}
