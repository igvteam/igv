/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature;

import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.List;

import static junit.framework.Assert.*;


/**
 * @author Jim Robinson
 * @date 5/10/12
 */
public class GFFTest{//} extends AbstractHeadlessTest{


    private List<Feature> getFeatures(String filePath) throws Exception{
        GFFParser parser = new GFFParser(filePath);
        BufferedReader br = new BufferedReader(new FileReader(filePath));
        List<Feature> features = parser.loadFeatures(br, null);
        br.close();
        return features;
    }
    /**
     * This test verifies that the attributes from column 9 are retained for a CDS feature that does not have a parent.
     * This bugfix was released with 2.1.12.
     *
     * @throws Exception
     */
    @Test
    public void testLoadSingleCDS() throws Exception {

        String gtfFile = TestUtils.DATA_DIR + "gtf/single_cds.gtf";
        List<Feature> features = getFeatures(gtfFile);

        assertEquals(1, features.size());
        BasicFeature bf = (BasicFeature) features.get(0);

        Exon exon = bf.getExons().get(0);
        assertNotNull(exon.getAttributes());
    }

    /**
     * Test that we parse the various spellings of "color" correctly
     * v2.2.0/2.2.1 had a bug where an attribute would be checked for
     * but a different one accessed
     * @throws Exception
     */
    @Test
    public void testGFFColorSpellings() throws Exception{
        String gffFile = TestUtils.DATA_DIR + "gff/color_sps.gff";
        List<Feature> features = getFeatures(gffFile);

        assertEquals("Error parsing certain features", 5, features.size());

        for(Feature f: features){
            BasicFeature bf = (BasicFeature) f;
            boolean haveColor = bf.getColor() != null;
            if(bf.hasExons()){
                for(Exon ex: bf.getExons()){
                    haveColor |= ex.getColor() != null;
                }
            }
            assertTrue(haveColor);
        }
    }


//    @Ignore
//    @Test
//    public void testLoadLargeFile() throws Exception{
//        //Very large file, as of this writing not in repo
//        String path = TestUtils.DATA_DIR + "gff/gene_coding_transcript_ncRNA_partiall_confirmed_plus.gff3";
//        GFFParser parser = new GFFParser(path);
//
//        BufferedReader br = new BufferedReader(new FileReader(path));
//        List<Feature> features = parser.loadFeatures(br, null);
//    }
}
