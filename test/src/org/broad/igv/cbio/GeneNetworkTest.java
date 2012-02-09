/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.cbio;

import org.apache.commons.collections.Predicate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.jgrapht.Graph;
import org.jgrapht.generate.WheelGraphGenerator;
import org.jgrapht.graph.Pseudograph;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import static junit.framework.Assert.*;

/**
 * Test our GeneNetwork class, which implements jgraphT graph interface.
 * {@see http://www.jgrapht.org/}.
 * <p/>
 * Notes:
 * Modification methods tend to return true for any modification,
 * and false for none. So graph.removeAll(nodeSet) would return
 * true if only half of the nodes were actually removed, although
 * odds are an exception would be thrown for any which could
 * not be removed.
 * <p/>
 * User: jacob
 * Date: 2012/02/02
 */
public class GeneNetworkTest {


    Genome genome;
    Graph<Vertex, BaseElement> graph;
    GraphMLExporter exporter;

    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();
        genome = TestUtils.loadGenome();
        graph = new Pseudograph<Vertex, BaseElement>(BaseElement.class);
        exporter = new GraphMLExporter();
    }

    @After
    public void tearDown() throws Exception {
        graph = null;
        exporter = null;
    }

    @Test
    public void testExportGraph() throws Exception {
        String filepath = TestUtils.DATA_DIR + "/out/test_cbio.graphml";
        String[] gene_list = new String[]{"EGFR", "BRCA1"};

        GeneNetwork geneNetwork = new GeneNetwork();
        geneNetwork.loadCBioLines(gene_list);

        exporter.exportGraph(filepath, geneNetwork, geneNetwork.factoryList);
    }

    @Test
    public void testLoadMutations() throws Exception {

        //String[] gene_list = new String[]{"TP53", "BRCA1", "KRAS"};
        String filepath = "test/data/out/test_cbio.graphml";
        String[] gene_list = new String[]{"EGFR", "BRCA1"};

        GeneNetwork geneNetwork = new GeneNetwork();
        //TODO We have defined loadCBioLines to return the number of vertices added
        //which is a bit arbitrary. may remove this later
        int initSize = geneNetwork.loadCBioLines(gene_list);

        //Load some mutation data
        String path = TestUtils.DATA_DIR + "/maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf.gz";
        ResourceLocator locator = new ResourceLocator(path);

        //IGV igv = TestUtils.startGUI();
        List<Track> tracks = new TrackLoader().load(locator, genome);

        for (Vertex v : geneNetwork.vertexSet()) {
            assertFalse(v.containsKey(GeneNetwork.PERCENT_MUTATED));
        }

        geneNetwork.collectMutationData(tracks);

        for (Vertex v : geneNetwork.vertexSet()) {
            assertTrue(v.containsKey(GeneNetwork.PERCENT_MUTATED));
        }

        exporter.exportGraph(filepath, geneNetwork, geneNetwork.factoryList);
    }

    @Test
    public void testFilterGraph() throws Exception {
        String filepath = "test/data/out/test_cbio.graphml";
        String[] gene_list = new String[]{"EGFR", "BRCA1"};
        final String key = GeneNetwork.COLUMN_NAMES[2];
        final String badval = "NA";

        GeneNetwork geneNetwork = new GeneNetwork();
        int initSize = geneNetwork.loadCBioLines(gene_list);

        Predicate has_evidence = new Predicate<BaseElement>() {
            public boolean evaluate(BaseElement object) {
                return object.containsKey(key)
                        && !badval.equals(object.get(key));
            }
        };


        assertTrue(geneNetwork.filterEdges(has_evidence));

        for (BaseElement e : geneNetwork.edgeSet()) {
            assertFalse(badval.equals(e.get(key)));
        }
    }


    /**
     * Test that graph is consistent (no duplicate vertices or edges).
     * Two vertices should be considered equal if their names are
     * equal. We give them different attributes
     *
     * @throws Exception
     */
    @Test
    public void testGraphVertices() throws Exception {
        KeyFactory factory = new KeyFactory("testGraphVertices");
        Random random = new Random(132294492231L);

        int numv = 20;
        String[] vertex_labels = new String[numv];
        String[] seed_attribs = new String[]{"attr1", "score", "squizzle", "thenadierwalzoftreachery"};
        String[] test_attribs = seed_attribs.clone();

        for (int ii = 0; ii < numv; ii++) {
            vertex_labels[ii] = "v" + ii;
            Vertex v = new Vertex(vertex_labels[ii], factory);
            assertTrue(graph.addVertex(v));

            for (String attr : seed_attribs) {
                v.put(attr, random.nextDouble());
            }
        }

        assertEquals(numv, graph.vertexSet().size());

        for (String label : vertex_labels) {
            Vertex v = new Vertex(label, factory);
            for (String attr : test_attribs) {
                v.put(attr, random.nextDouble());
            }
            assertFalse(graph.addVertex(v));
        }
        assertEquals(numv, graph.vertexSet().size());
    }

    @Test
    public void testPruneGraph() throws Exception {
        int size = 10;
        int exp_edges = size - 1 + size - 1;
        KeyFactory nodeKF = new KeyFactory("testPruneGraph");
        KeyFactory edgeKF = new KeyFactory("testPruneGraphEdge");
        GeneNetwork.SimpleVertexFactory vFactory = new GeneNetwork.SimpleVertexFactory(nodeKF);
        GeneNetwork.SimpleEdgeFactory eFactory = new GeneNetwork.SimpleEdgeFactory(edgeKF);
        GeneNetwork graph = new GeneNetwork(eFactory);


        //Bicycle wheel. Hub node is connected to all elements, outer elements
        //conected in ring.
        WheelGraphGenerator<Vertex, BaseElement> wgg = new WheelGraphGenerator<Vertex, BaseElement>(size);
        wgg.generateGraph(graph, vFactory, null);
        assertEquals(size, graph.vertexSet().size());
        assertEquals(exp_edges, graph.edgeSet().size());

        //Find first spoke vertex, isolate it.
        Vertex toRem = null;
        for (Vertex v : graph.vertexSet()) {
            if (graph.edgesOf(v).size() < size && graph.edgesOf(v).size() == 3) {
                //Spoke vertex
                toRem = v;
                break;
            }
        }

        Set<BaseElement> edges = new HashSet<BaseElement>();
        for (BaseElement e : graph.edgesOf(toRem)) {
            edges.add(e);
        }
        assertTrue(graph.removeAllEdges(edges));

        //At this point, toRem should be isolated
        assertEquals(0, graph.edgesOf(toRem).size());
        assertTrue(graph.pruneGraph());

        assertEquals(size - 1, graph.vertexSet().size());
        assertEquals(exp_edges - 3, graph.edgeSet().size());
        assertFalse(graph.containsVertex(toRem));

        int trials = 10;
        for (int ii = 0; ii < trials; ii++) {
            assertFalse(graph.pruneGraph());
        }

    }


    /*
    for(Vertex v: rejected){
        assertTrue(geneNetwork.removeVertex(v));
    }
    //assertTrue(rejected.size() > 0);
    assertEquals(initSize - rejected.size(), geneNetwork.vertexSet().size());



    List<Track> tracks = IGV.getInstance().getAllTracks(false);
    for(Track track: tracks){
        if (track.getTrackType() == TrackType.MUTATION) {
            track.getRegionScore()?
        }
    }


    MutationParser parser = new MutationParser();
    Map<String, List<org.broad.tribble.Feature>> mutations =  parser.loadMutations(locator, genome);
    System.out.println(mutations.size());

    for(String s: mutations.keySet()){
        List<Feature> features = mutations.get(s);
        for(Feature f: features){
            System.out.println(f.getClass());
            //This doesn't work for obvious reasons, but the flow should be clear.
            //Load the list of features from FeatureDB
            //FeatureUtils.getFeatureAt(f.getStart(), f.getEnd() - f.getStart(), graph.vertexSet());
        }

    }
    */

}
