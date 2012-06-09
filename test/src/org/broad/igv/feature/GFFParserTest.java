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

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;


/**
 * @author Jim Robinson
 * @date 5/10/12
 */
public class GFFParserTest {

    /**
     * This test verifies that the attributes from column 9 are retained for a CDS feature that does not have a parent.
     * This bugfix was released with 2.1.12.
     *
     * @throws Exception
     */
    @Test
    public void testLoadSingleCDS() throws Exception {

        String gtfFile = TestUtils.DATA_DIR + "gtf/single_cds.gtf";
        GFFParser parser = new GFFParser(gtfFile);
        Genome genome = null;
        BufferedReader br = new BufferedReader(new FileReader(gtfFile));
        List<Feature> features = parser.loadFeatures(br, genome);
        br.close();

        assertEquals(1, features.size());
        BasicFeature bf = (BasicFeature) features.get(0);

        Exon exon = bf.getExons().get(0);
        assertNotNull(exon.getAttributes());


    }
}
