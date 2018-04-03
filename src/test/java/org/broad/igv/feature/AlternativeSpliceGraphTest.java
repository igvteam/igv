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

package org.broad.igv.feature;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.jgrapht.ext.GmlExporter;
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
public class AlternativeSpliceGraphTest extends AbstractHeadlessTest {

    @Test
    public void initTest() throws Exception {

        String nm = "CLIP1";
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
