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

package org.broad.igv.cbio;

import biz.source_code.base64Coder.Base64Coder;
import com.google.common.base.Predicate;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.*;
import org.broad.igv.util.collections.CollUtils;
import org.jgrapht.EdgeFactory;
import org.jgrapht.VertexFactory;
import org.jgrapht.generate.WheelGraphGenerator;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import javax.xml.parsers.DocumentBuilderFactory;
import java.io.*;
import java.util.Collection;
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
public class GeneNetworkTest extends AbstractHeadlessTest {

    private static String testpath = TestUtils.DATA_DIR + "xml/tp53network.xml";
    private static String drugTestPath = TestUtils.DATA_DIR + "xml/EGFR_withdrugs.xml.gz";
    private GeneNetwork network;

    @Before
    public void setUp() throws Exception {
        network = new GeneNetwork();
    }

    @After
    public void tearDown() throws Exception {
        network = null;
        GeneNetwork.BASE_URL = GeneNetwork.REAL_URL;
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

    /**
     * Test that we don't filter out drug nodes
     *
     * @throws Exception
     */
    @Test
    public void testFilterWithDrugs() throws Exception {
        assertTrue("Failed to load network", network.loadNetwork(drugTestPath) > 0);
        Collection<Node> nonGenesBeforeFilter = CollUtils.filter(network.vertexSet(), GeneNetwork.isNotGene);
        assertTrue("Bad test setup, no non-gene nodes", nonGenesBeforeFilter.size() > 0);

        doTestAnnotation(network);
        int genesRemoved = network.filterGenesRange(GeneNetwork.PERCENT_MUTATED, 0, 10.0f);
        assertTrue("Bad test setup, Filter didn't remove any genes", genesRemoved > 0);

        Set<Node> nonGenesAfterFilter =
                new HashSet<Node>(CollUtils.filter(network.vertexSet(), GeneNetwork.isNotGene));
        for (Node nonGene : nonGenesBeforeFilter) {
            String msg = String.format("Filtered out node %s which we shouldn't have", GeneNetwork.getNodeKeyData(nonGene, GeneNetwork.LABEL));
            assertTrue(msg, nonGenesAfterFilter.contains(nonGene));
        }

    }

    private static Predicate evenPred = new Predicate() {

        public boolean apply(Object object) {
            Node node = (Node) object;
            NamedNodeMap map = node.getAttributes();
            if (map == null) {
                return false;
            }
            int id = Integer.parseInt(map.getNamedItem("id").getTextContent());
            return id % 2 == 0;
        }
    };

    @Test
    public void testReset() throws Exception {
        network.loadNetwork(testpath);
        int initSize = network.vertexSet().size();
        int nonQuery = network.filterGenes(GeneNetwork.inQuery);
        assertTrue("Bad test setup, all genes in query", nonQuery > 0);
        assertEquals(initSize - network.vertexSet().size(), nonQuery);

        network.reset();
        assertEquals("Resetting network did not bring back everything", initSize, network.vertexSet().size());
    }

    @Test
    public void testFilterKeepsQuery() throws Exception {
        network.loadNetwork(testpath);

        network.filterGenes(GeneNetwork.inQuery);

        Collection<Node> queryGenes = network.geneVertexes();
        network.reset();

        int numRemoved = network.filterGenes(evenPred);
        assertTrue(numRemoved > 0);
        Collection<Node> randomGenes = network.geneVertexes();
        Set<Node> randomGeneSet = new HashSet<Node>(randomGenes);

        for (Node queryGene : queryGenes) {
            assertTrue(randomGeneSet.contains(queryGene));
        }
    }

    @Test
    public void testFilterKeepsEdges() throws Exception {
        network.loadNetwork(testpath);
        int initSize = network.vertexSet().size();

        int numRemoved = network.filterGenes(evenPred);
        assertTrue(numRemoved > 0);

        //Test that we can get the filtered edges of a node
        Set<Node> keptNodes = new HashSet<Node>();
        for (Node n : network.geneVertexes()) {
            for (Node e : network.edgesOf(n)) {
                keptNodes.add(network.getEdgeSource(e));
                keptNodes.add(network.getEdgeTarget(e));
            }
        }
        assertEquals(network.geneVertexes().size(), keptNodes.size());
        assertTrue("Filtering not performed", keptNodes.size() < initSize);
    }


    @Test
    public void testFilterOnEvidence() throws Exception {
        GeneNetwork geneNetwork = new GeneNetwork();
        assertTrue(geneNetwork.loadNetwork(TestUtils.DATA_DIR + "xml/egfr_brca1.xml.gz") > 0);

        final String badname = "NA";

        Predicate<Node> has_evidence = new Predicate<Node>() {
            public boolean apply(Node object) {
                String label = GeneNetwork.getNodeKeyData(object, "EXPERIMENTAL_TYPE");
                return label != null && !label.equals(badname);
            }
        };

        //Perform the filtering.
        //This is true if any modifications are made.
        assertTrue(geneNetwork.filterEdges(has_evidence) > 0);

        for (Node e : geneNetwork.edgeSet()) {
            assertTrue(has_evidence.apply(e));
        }
    }

    @Test
    public void testOutputNoGzip() throws Exception {
        String networkPath = TestUtils.DATA_DIR + "xml/egfr_brca1.xml.gz";
        //String networkPath = testpath;
        assertTrue(network.loadNetwork(networkPath) > 0);
        String outPath = TestUtils.DATA_DIR + "out/test.xml";
        tstOutputNetwork(network, outPath);
    }

    @Test
    public void testOutputGzip() throws Exception {
        String networkPath = TestUtils.DATA_DIR + "xml/egfr_brca1.xml.gz";
        assertTrue(network.loadNetwork(networkPath) > 0);
        String outPath = TestUtils.DATA_DIR + "out/test.xml.gz";
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

    private void doTestAnnotation(GeneNetwork network) throws Exception {

        //Load some tracks
        String dataPath = TestUtils.DATA_DIR + "seg/Broad.080528.subtypes.seg.gz";
        ResourceLocator locator = new ResourceLocator(dataPath);
        List<Track> tracks = new TrackLoader().load(locator, genome);
        network.annotateAll(tracks);
    }


    @Test
    public void testAnnotateAll() throws Exception {

        String networkPath = TestUtils.DATA_DIR + "xml/egfr_brca1.xml.gz";
        assertTrue(network.loadNetwork(networkPath) > 0);
        doTestAnnotation(network);

        //Check data
        Set<Node> nodes = network.vertexSet();
        for (Node node : nodes) {
            for (String key : GeneNetwork.attributeMap.keySet()) {
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
        for (String key : GeneNetwork.attributeMap.keySet()) {
            String data = GeneNetwork.getNodeAttrValue(gml, "id", key);
            assertNotNull(data);
        }

    }

    @Test
    public void testPruneGraph() throws Exception {
        int size = 10;
        int exp_edges = size - 1 + size - 1;
        Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        SimpleVertexFactory vFactory = new SimpleVertexFactory(document);
        SimpleEdgeFactory eFactory = new SimpleEdgeFactory(document);
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

        assertEquals(size - 1, graph.vertexSet().size());
        assertEquals(exp_edges - 3, graph.edgeSet().size());

        assertFalse(graph.containsVertex(toRem));
    }

    @Test
    public void testOutputForcBioView() throws Exception {
        assertTrue(network.loadNetwork(testpath) > 0);
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

    /**
     * Generates nodes with unique IDs, all with <data key="type">Protein</data>
     * attributes
     */
    public static class SimpleVertexFactory implements VertexFactory<Node> {
        private long vCounter = 0;

        private final Document document;

        public SimpleVertexFactory(Document document) {
            this.document = document;
        }

        public Node createVertex() {
            Element v = document.createElement("v" + vCounter);
            v.setAttribute("id", "" + vCounter);

            //<data key="type">Protein</data>
            Element type = document.createElement("data");
            type.setAttribute("key", "type");
            type.setTextContent("Protein");
            v.appendChild(type);
            vCounter++;
            return v;
        }
    }

    public static class SimpleEdgeFactory implements EdgeFactory<Node, Node> {
        private long eCounter = 0;
        private final Document document;

        public SimpleEdgeFactory(Document document) {
            this.document = document;
        }

        public Node createEdge(Node sourceVertex, Node targetVertex) {
            Element e = document.createElement("e" + eCounter);
            e.setAttribute("source", sourceVertex.getAttributes().getNamedItem("id").toString());
            e.setAttribute("target", targetVertex.getAttributes().getNamedItem("id").toString());
            eCounter++;
            return e;
        }
    }

}
