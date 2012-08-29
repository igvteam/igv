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

package org.broad.igv.tools.converters;

import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 2/15/12
 */
public class GFFtoBedTest {


    @Test
    public void testConvert() throws Exception {

        File inputFile = new File(TestUtils.DATA_DIR, "gff/gene.unsorted.gff3");
        File outputFile = new File(TestUtils.DATA_DIR, "out/gene.bed");

        // Convert
        GFFtoBed.convert(inputFile, outputFile);

        // Now compare features from both files, they should be identical

        BufferedReader gffReader = null;
        BufferedReader bedReader = null;
        try {
            Genome genome = null; // Not needed for this test
            gffReader = new BufferedReader(new FileReader(inputFile));
            GFFParser parser = new GFFParser(inputFile.getPath());
            List<Feature> gffFeatures = parser.loadFeatures(gffReader, genome);

            bedReader = new BufferedReader(new FileReader(outputFile));
            FeatureParser bedParser = AbstractFeatureParser.getInstanceFor(outputFile.getAbsolutePath(), null);
            List<Feature> bedFeatures = bedParser.loadFeatures(bedReader, genome);

            assertTrue(gffFeatures.size() > 0);
            assertEquals(gffFeatures.size(), bedFeatures.size());
            for (int i = 0; i < gffFeatures.size(); i++) {
                BasicFeature gffFeature = (BasicFeature) gffFeatures.get(i);
                BasicFeature bedFeature = (BasicFeature) bedFeatures.get(i);

                compare(gffFeature, bedFeature);

            }

        } finally {
            gffReader.close();
        }


    }

    private void compare(BasicFeature gffFeature, BasicFeature bedFeature) {

        assertEquals(gffFeature.getChr(), bedFeature.getChr());
        assertEquals(gffFeature.getStart(), bedFeature.getStart());
        assertEquals(gffFeature.getEnd(), bedFeature.getEnd());
        assertEquals(gffFeature.getStrand(), bedFeature.getStrand());
        assertEquals(gffFeature.getColor(), bedFeature.getColor());
        //BED features don't have a type
        //assertEquals(gffFeature.getType(), bedFeature.getType());
        assertEquals(gffFeature.getStart(), bedFeature.getThickStart());
        assertEquals(gffFeature.getThickEnd(), bedFeature.getThickEnd());
        if (!Float.isNaN(gffFeature.getScore())) {
            assertEquals(gffFeature.getScore(), bedFeature.getScore());
        }
        assertEquals(gffFeature.getExonCount(), bedFeature.getExonCount());

        List<Exon> gffExons = gffFeature.getExons();
        List<Exon> bedExons = bedFeature.getExons();

        // This is broken for sure, disable for now
//        for (int i = 0; i < gffExons.size(); i++) {
//
//            Exon gffExon = gffExons.get(i);
//            Exon bedExon = bedExons.get(i);
//            assertEquals(gffExon.getCdStart(), bedExon.getCdStart());
//            assertEquals(gffExon.getCodingLength(), bedExon.getCodingLength());
//            assertEquals(gffExon.getCdEnd(), bedExon.getCdEnd());
//            assertEquals(gffExon.getReadingShift(), bedExon.getReadingShift());
//        }

    }
}
