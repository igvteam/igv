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

package org.broad.igv.feature;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.jgrapht.ext.GmlExporter;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/03/28
 */
@Ignore
public class AlternativeSpliceGraphTest {
    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();
        Genome genome = TestUtils.loadGenome();
    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void initTest() throws Exception {

        String nm = "EGFR";
        List<NamedFeature> n_features = FeatureDB.getFeaturesList(nm, Integer.MAX_VALUE, false);
        List<IGVFeature> features = new ArrayList<IGVFeature>(n_features.size());
        int maxNumExons = -1;
        for (NamedFeature feat : n_features) {
            if (feat.getName().equals(nm) && feat instanceof IGVFeature) {
                IGVFeature iFeat = (IGVFeature) feat;
                if (iFeat.getExons().size() > maxNumExons) {
                    maxNumExons = iFeat.getExons().size();
                }
                features.add(iFeat);
            }
        }

        AlternativeSpliceGraph graph = new AlternativeSpliceGraph(features);
        assertTrue(graph.vertexSet().size() >= maxNumExons);

        GmlExporter exporter = new GmlExporter();
        String outfile = TestUtils.DATA_DIR + "testas.gml";
        FileWriter writer = new FileWriter(outfile);
        exporter.export(writer, graph);


    }

}
