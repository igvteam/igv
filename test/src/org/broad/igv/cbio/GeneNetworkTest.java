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

import biz.source_code.base64Coder.Base64Coder;
import org.apache.commons.collections.Predicate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.*;
import org.jgrapht.EdgeFactory;
import org.jgrapht.VertexFactory;
import org.jgrapht.generate.WheelGraphGenerator;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import javax.imageio.metadata.IIOMetadataNode;
import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

/**
 * Test our GeneNetwork class, which implements jgraphT graph interface.
 * See  http://www.jgrapht.org.
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

    private static String testpath = TestUtils.DATA_DIR + "/tp53network.xml";
    private GeneNetwork network;

    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();
        network = new GeneNetwork();
    }

    @After
    public void tearDown() throws Exception {
        network = null;
        TestUtils.clearOutputDir();
    }

    public static void main(String[] args) throws IOException {
        String netpath = testpath;
        if (args != null && args.length > 0) {
            netpath = args[0];
        }

        GeneNetwork network = new GeneNetwork();
        network.loadNetwork(netpath);
        String viewPath = network.outputForcBioView();
        BrowserLauncher.openURL("file://" + viewPath);
    }

    @Test
    public void testLoadLocal() throws Exception {
        assertTrue("Failed to load network", network.loadNetwork(testpath) > 0);
    }

    @Test
    public void testFilter() throws Exception {
        Predicate tPred = new Predicate() {

            public boolean evaluate(Object object) {
                Node node = (Node) object;
                NamedNodeMap map = node.getAttributes();
                if (map == null) {
                    return false;
                }
                int id = Integer.parseInt(map.getNamedItem("id").getTextContent());
                return id % 2 == 0;
            }
        };

        network.loadNetwork(testpath);
        boolean removed = network.filterNodes(tPred) > 0;
        assertTrue(removed);
    }

    /**
     * Load some data from cbio.
     * Checks that we are looking at the right urls
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void testDownloadCBIO() throws Exception {
        String[] gene_list = new String[]{"egfr", "brca1", "jun"};
        GeneNetwork anno = GeneNetwork.getFromCBIO(Arrays.asList(gene_list));
        assertNotNull(anno);
    }

    @Test
    public void testOutputNoGzip() throws Exception {
        String networkPath = TestUtils.DATA_DIR + "/egfr_brca1.xml.gz";
        //String networkPath = testpath;
        assertTrue(network.loadNetwork(networkPath) > 0);
        String outPath = TestUtils.DATA_DIR + "/out/test.xml";
        tstOutputNetwork(network, outPath);
    }

    @Test
    public void testOutputGzip() throws Exception {
        String networkPath = TestUtils.DATA_DIR + "/egfr_brca1.xml.gz";
        assertTrue(network.loadNetwork(networkPath) > 0);
        String outPath = TestUtils.DATA_DIR + "/out/test.xml.gz";
        tstOutputNetwork(network, outPath);
    }


    public void tstOutputNetwork(GeneNetwork network, String outPath) throws Exception {
        Set<Node> nodes = network.vertexSet();
        Set<String> nodeNames = new HashSet<String>();

        for (Node node : nodes) {
            nodeNames.add(GeneNetwork.getNodeKeyData(node, "label"));
        }

        assertTrue(network.exportGraph(outPath) > 0);


        GeneNetwork at = new GeneNetwork();
        assertTrue(at.loadNetwork(outPath) > 0);

        //Check that node set matches
        Set<Node> outNodes = at.vertexSet();
        assertEquals("Output has a different number of nodes than input", nodes.size(), outNodes.size());
        for (Node oNode : outNodes) {
            String nodeName = GeneNetwork.getNodeKeyData(oNode, "label");
            assertTrue(nodeNames.contains(nodeName));
        }

    }


    @Test
    public void testAnnotateAll() throws Exception {
        TestUtils.setUpHeadless();
        Genome genome = TestUtils.loadGenome();

        String networkPath = TestUtils.DATA_DIR + "/egfr_brca1.xml.gz";
        assertTrue(network.loadNetwork(networkPath) > 0);

        //Load some tracks
        String dataPath = TestUtils.DATA_DIR + "/seg/Broad.080528.subtypes.seg.gz";
        ResourceLocator locator = new ResourceLocator(dataPath);
        List<Track> tracks = new TrackLoader().load(locator, genome);
        network.annotateAll(tracks);

        //Check data
        Set<Node> nodes = network.vertexSet();
        for (Node node : nodes) {
            for (String key : GeneNetwork.attribute_map.keySet()) {
                String data = GeneNetwork.getNodeKeyData(node, key);
                String name = GeneNetwork.getNodeKeyData(node, GeneNetwork.LABEL);
                if (!"CHMP3".equalsIgnoreCase(name)) {
                    assertNotNull(data);
                }
            }
        }

        //Check schema
        Document doc = network.createDocument();
        Node gml = doc.getFirstChild();
        for (String key : GeneNetwork.attribute_map.keySet()) {
            String data = GeneNetwork.getNodeAttrValue(gml, "id", key);
            assertNotNull(data);
        }

    }

    @Test
    public void testFilterGraph() throws Exception {
        GeneNetwork geneNetwork = new GeneNetwork();
        assertTrue(geneNetwork.loadNetwork(TestUtils.DATA_DIR + "/egfr_brca1.xml.gz") > 0);

        final String badname = "NA";

        Predicate<Node> has_evidence = new Predicate<Node>() {
            public boolean evaluate(Node object) {
                String label = GeneNetwork.getNodeKeyData(object, "EXPERIMENTAL_TYPE");
                return label != null && !label.equals(badname);
            }
        };

        //Perform the filtering.
        //This is true if any modifications are made.
        assertTrue(geneNetwork.filterEdges(has_evidence) > 0);

        geneNetwork.finalizeFilters();
        for (Node e : geneNetwork.edgeSet()) {
            assertTrue(has_evidence.evaluate(e));
        }
    }


    @Test
    public void testPruneGraph() throws Exception {
        int size = 10;
        int exp_edges = size - 1 + size - 1;
        SimpleVertexFactory vFactory = new SimpleVertexFactory();
        SimpleEdgeFactory eFactory = new SimpleEdgeFactory();
        GeneNetwork graph = new GeneNetwork(eFactory);


        //Bicycle wheel. Hub node is connected to all elements, outer elements
        //conected in ring.
        WheelGraphGenerator<Node, Node> wgg = new WheelGraphGenerator<Node, Node>(size);
        wgg.generateGraph(graph, vFactory, null);
        assertEquals(size, graph.vertexSet().size());
        assertEquals(exp_edges, graph.edgeSet().size());

        //Find first spoke vertex, isolate it.
        Node toRem = null;
        for (Node v : graph.vertexSet()) {
            if (graph.edgesOf(v).size() < size && graph.edgesOf(v).size() == 3) {
                //Spoke vertex
                toRem = v;
                break;
            }
        }

        Set<Node> edges = new HashSet<Node>();
        for (Node e : graph.edgesOf(toRem)) {
            edges.add(e);
        }
        assertTrue(graph.removeAllEdges(edges));

        //At this point, toRem should be isolated
        assertEquals(0, graph.edgesOf(toRem).size());
        assertTrue(graph.pruneGraph());
        graph.finalizeFilters();

        assertEquals(size - 1, graph.vertexSet().size());
        assertEquals(exp_edges - 3, graph.edgeSet().size());
        assertFalse(graph.containsVertex(toRem));
    }

    @Test
    public void testOutputForcBioView() throws Exception {
        assertTrue(network.loadNetwork(TestUtils.DATA_DIR + "/tp53network.xml") > 0);
        String outPath = network.outputForcBioView();

        //Now attempt to read back in
        //Note that this will not work if the output is in plaintext,
        //because it contains the XML header in the middle of the document
        Document inDoc = Utilities.createDOMDocumentFromXmlStream(new FileInputStream(outPath));
        String b64data = inDoc.getElementsByTagName("textarea").item(0).getTextContent().trim();
        byte[] gzippedInput = Base64Coder.decode(b64data);

        BufferedReader bufIn = null;
        try {
            InputStream plainData = new GZIPInputStream(new ByteArrayInputStream(gzippedInput));
            bufIn = new BufferedReader(new InputStreamReader(plainData));
        } catch (IOException e) {
            bufIn = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(gzippedInput)));
        }

        int count = 0;
        String[] outLines = Utilities.getString(network.createDocument()).split(FileUtils.LINE_SEPARATOR);
        String line;
        while ((line = bufIn.readLine()) != null) {
            assertEquals(outLines[count], line);
            count++;
        }
    }

    public static class SimpleVertexFactory implements VertexFactory<Node> {
        private long vCounter = 0;

        public SimpleVertexFactory() {
        }

        public Node createVertex() {
            Element v = new IIOMetadataNode("v" + vCounter);
            v.setAttribute("id", "" + vCounter);
            vCounter++;
            return v;
        }
    }

    public static class SimpleEdgeFactory implements EdgeFactory<Node, Node> {
        private long eCounter = 0;

        public SimpleEdgeFactory() {
        }

        public Node createEdge(Node sourceVertex, Node targetVertex) {
            Element e = new IIOMetadataNode("e" + eCounter);
            e.setAttribute("source", sourceVertex.getAttributes().getNamedItem("id").toString());
            e.setAttribute("target", targetVertex.getAttributes().getNamedItem("id").toString());
            eCounter++;
            return e;
        }
    }


}
