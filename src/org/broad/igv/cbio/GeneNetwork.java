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
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.*;
import org.jgrapht.EdgeFactory;
import org.jgrapht.graph.Pseudograph;
import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * Fetch a gene network from cBio portal. WORK IN PROGRESS.
 * User: jacob
 * Date: 2012/01/31
 */
public class GeneNetwork extends Pseudograph<Node, Node> {

    private Logger logger = Logger.getLogger(GeneNetwork.class);

    public static final String NODE_TAG = "node";
    public static final String EDGE_TAG = "edge";
    public static final String KEY = "key";
    public static final String LABEL = "label";


    /**
     * URL that cbio will use when service is released
     */
    public static String REAL_URL = "http://www.cbioportal.org/public-portal/webservice.do";
    /**
     * URL they use for testing
     */
    public static String TEST_URL = "http://awabi.cbio.mskcc.org/public-portal/network.do";
    public static String BASE_URL = TEST_URL;
    //    static{
//        InputStream is = NetworkAnnotator.class.getResourceAsStream("resources/url.txt");
//        BufferedReader br = new BufferedReader(new InputStreamReader(is));
//        try{
//            BASE_URL = br.readLine();
//        }catch(IOException e){
//            logger.error("url resource not found, defaulting to " + TEST_URL);
//            BASE_URL = TEST_URL;
//        }
//    }
    private static final String common_parms = "format=gml&gzip=on";
    private static final String GENE_LIST = "gene_list";

    private Map<String, Node> nodeTable;

    private List<Node> schema = new ArrayList<Node>();
    private NamedNodeMap graphAttr;
    private Document origDocument;


    static Map<String, RegionScoreType> attribute_map = new HashMap();

    static {
        attribute_map.put("PERCENT_MUTATED", RegionScoreType.MUTATION_COUNT);
        attribute_map.put("PERCENT_CNA_AMPLIFIED", RegionScoreType.AMPLIFICATION);
        attribute_map.put("PERCENT_CNA_HOMOZYGOUSLY_DELETED", RegionScoreType.DELETION);
    }

    public static final String PERCENT_ALTERED = "PERCENT_ALTERED";

//    protected KeyFactory edgeKeyFactory;
//    protected KeyFactory vertexKeyFactory;
//    protected List<KeyFactory> factoryList;

    public GeneNetwork() {
        this(Node.class);
    }

    public GeneNetwork(EdgeFactory edgeFactory) {
        super(edgeFactory);
    }

    public GeneNetwork(Class clazz) {
        super(clazz);
//
//        edgeKeyFactory = new KeyFactory("edge");
//        vertexKeyFactory = new KeyFactory("node");
//        factoryList = Arrays.asList(edgeKeyFactory, vertexKeyFactory);
    }

//    public List<KeyFactory> getFactoryList() {
//        return this.factoryList;
//    }

    public boolean filterNodes(Predicate predicate) {
        Set<Node> rejected = new HashSet<Node>(this.vertexSet().size());
        for (Node v : this.vertexSet()) {
            if (!predicate.evaluate(v)) {
                rejected.add(v);
            }
        }
        return this.removeAllVertices(rejected);
    }

    public boolean filterNodesRange(final String key, final float min, final float max) {
        Predicate<Node> pred = new Predicate<Node>() {
            @Override
            public boolean evaluate(Node object) {
                String sval = getNodeKeyData(object, key);
                if (sval != null) {
                    float fval = Float.parseFloat(sval);
                    return fval >= min && fval <= max;
                }
                return false;
            }
        };

        return this.filterNodes(pred);
    }

    public boolean filterEdges(Predicate predicate) {
        Set<Node> rejected = new HashSet<Node>(this.edgeSet().size());
        for (Node e : this.edgeSet()) {
            if (!predicate.evaluate(e)) {
                rejected.add(e);
            }
        }
        return this.removeAllEdges(rejected);
    }

    public boolean pruneGraph() {
        Predicate min_connections = new Predicate<Node>() {
            public boolean evaluate(Node object) {
                return edgesOf(object).size() >= 1;
            }
        };
        return this.filterNodes(min_connections);
    }

    public static GeneNetwork getFromCBIO(Iterable<String> geneList) {
        String query = HttpUtils.buildURLString(geneList, "+");
        String url = BASE_URL + "?" + GENE_LIST + "=" + query + "&" + common_parms;
        GeneNetwork network = new GeneNetwork();
        if (network.loadNetwork(url) > 0) {
            return network;
        } else {
            return null;
        }
    }

    /**
     * Load graphml data from provide path.
     *
     * @param path
     */
    public int loadNetwork(String path) {
        String error = null;
        int numNodes = -1;
        try {
            InputStream cbioStream = ParsingUtils.openInputStreamGZ(new ResourceLocator(path));
            Document document = Utilities.createDOMDocumentFromXmlStream(cbioStream);
            this.origDocument = document;

            //Read schema from top and save it
            addToSchema(document.getElementsByTagName("key"));

            graphAttr = document.getElementsByTagName("graph").item(0).getAttributes();


            NodeList nodes = document.getElementsByTagName(NODE_TAG);
            //Generate the graph itself. First add the nodes, then the edges
            nodeTable = new HashMap<String, Node>(nodes.getLength());
            for (int nn = 0; nn < nodes.getLength(); nn++) {
                Node node = nodes.item(nn);
                String label = node.getAttributes().getNamedItem("id").getTextContent();
                nodeTable.put(label, node);
                this.addVertex(node);
            }

            NodeList edges = document.getElementsByTagName(EDGE_TAG);
            for (int ee = 0; ee < edges.getLength(); ee++) {
                Node edge = edges.item(ee);
                NamedNodeMap attrs = edge.getAttributes();
                String source = attrs.getNamedItem("source").getTextContent();
                String target = attrs.getNamedItem("target").getTextContent();
                this.addEdge(nodeTable.get(source), nodeTable.get(target), edge);
            }

            numNodes = this.vertexSet().size();

        } catch (IOException e) {
            error = e.getMessage();
        } catch (ParserConfigurationException e) {
            error = e.getMessage();
        } catch (SAXException e) {
            error = e.getMessage();
        } finally {
            if (error != null) {
                logger.error(error);
                return -1;
            } else {
                return numNodes;
            }
        }
    }

    /**
     * Add schema information for the provided datakeys.
     * They will all be set as the provided dataType (string, double, float, etc.)
     * and graph element
     *
     * @param dataKeys
     * @param dataType Legal values are long, integer, float, double, boolean, string. Case insensitive.
     *                 All dataKeys will be set to the provided type.
     * @param typeFor  Legal values are "all", "graph", "node", "edge"
     */
    private void addSchema(Collection<String> dataKeys, String dataType, String typeFor) {
        Element key;

        for (String dataKey : dataKeys) {
            key = this.origDocument.createElement("key");
            key.setAttribute("id", dataKey);
            //TODO id is supposed to unique, attr.name human readable.
            //Not quite sure of any reason they can't be the same.
            key.setAttribute("attr.name", dataKey);
            key.setAttribute("attr.type", dataType.toLowerCase());

            if (typeFor != null) {
                key.setAttribute("for", typeFor);
            }
            schema.add(key);
        }
    }


    private void addToSchema(NodeList keys) {
        for (int kk = 0; kk < keys.getLength(); kk++) {
            schema.add(keys.item(kk));
        }
    }

    /**
     * The the value of a child node by the key.
     * If there are multiple matches, the first is returned.
     * Search is not recursive.
     * <p/>
     * <p/>
     * Example: Say that node has the following XML
     * "&lt;node id="3725"/&gt;
     * &lt;data key="label"&gt;JUN&lt;/data&gt;
     * &lt;data key="type"&gt;Protein&lt;/data&gt;
     * &lt;data key="RELATIONSHIP_XREF"&gt;HGNC:JUN;Entrez Gene:3725&lt;/data&gt;
     * &lt;data key="IN_QUERY"&gt;false&lt;/data&gt;
     * &lt;/node&gt;"
     * So getNodeKeyData(node, "key", "label") returns "JUN".
     *
     * @param node
     * @param attrName
     * @param attrValue
     * @return String value of key found. null if not found
     */
    public static String getNodeAttrValue(Node node, String attrName, String attrValue) {
        NodeList elements = node.getChildNodes();
        for (int ee = 0; ee < elements.getLength(); ee++) {
            Node el = elements.item(ee);
            try {
                NamedNodeMap map = el.getAttributes();
                Node label = map.getNamedItem(attrName);
                String textContent = label.getTextContent();
                if (textContent.compareToIgnoreCase(attrValue) == 0) {
                    return el.getTextContent();
                }
            } catch (NullPointerException e) {
                //In general these get hit due to newlines and such
                //We simply skip
                continue;
            }
        }
        return null;
    }

    /**
     * Gets the value of a child node with "key" attribute
     * equal to {@code key} parameter.
     * Equal to getNodeAttrValue(node, NetworkAnnotator.KEY, key);
     *
     * @param node Node to search
     * @param key  Key to search for
     * @return String value of key found. null if not found
     * @see #getNodeAttrValue
     */
    public static String getNodeKeyData(Node node, String key) {
        return getNodeAttrValue(node, KEY, key);
    }

    public Document createDocument() {
        try {
            // Create a DOM document
            DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document document = documentBuilder.newDocument();
            document.setStrictErrorChecking(false);

            // Global root element
            Element globalElement = document.createElement("graphml");

            //Add schema and attributes
            for (Node node : schema) {
                globalElement.appendChild(node);
            }

            Element graphEl = document.createElement("graph");

            for (int aa = 0; aa < graphAttr.getLength(); aa++) {
                Node attr = graphAttr.item(aa);
                graphEl.setAttribute(attr.getNodeName(), attr.getTextContent());
            }
            //Add nodes and edges
            for (Node v : this.vertexSet()) {
                graphEl.appendChild(v);
            }

            for (Node e : this.edgeSet()) {
                graphEl.appendChild(e);
            }

            globalElement.appendChild(graphEl);
            document.appendChild(globalElement);

            return document;
        } catch (Exception e) {
            throw new RuntimeException("Error outputting graph");
        }
    }

    public static int writeEncodedString(String string, OutputStream outputStream, boolean gzip, boolean base64encode) throws IOException {

        byte[] byteData;

        if (gzip) {
            ByteArrayOutputStreamChild gmlByteStream = new ByteArrayOutputStreamChild(string.length() / 20);
            GZIPOutputStream gzipOutputStream = new GZIPOutputStream(gmlByteStream);
            gzipOutputStream.write(string.getBytes());
            gzipOutputStream.finish();
            byteData = gmlByteStream.getBuf();
            gmlByteStream.close();
        } else {
            byteData = string.getBytes();
        }


        int count = 0;
        if (base64encode) {
            char[] gmlData = Base64Coder.encode(byteData);
            for (char c : gmlData) {
                outputStream.write(c);
                count++;
            }
        } else {
            outputStream.write(byteData);
            outputStream.flush();
            count += byteData.length;
        }
        outputStream.flush();
        return count;
    }


    /**
     * Write document to XML at outputFile. File is deleted if there
     * is an error writing out. If the outputFile has a .gz extension,
     * the output is gzipped.
     *
     * @param outputFile
     * @return success
     * @throws java.io.IOException
     */
    public int exportGraph(String outputFile) throws IOException {
        boolean gzip = outputFile.endsWith(".gz");

        String xmlString = Utilities.getString(this.createDocument());

        OutputStream stream = new FileOutputStream(outputFile);
        int count = writeEncodedString(xmlString, stream, gzip, false);

        stream.flush();
        stream.close();

        return count;
    }

    public String outputForcBioView() throws IOException {
        String outPath = null;
        BufferedReader reader = null;
        OutputStream outputStream = null;
        try {
            File temp = File.createTempFile("cbio", ".html");
            temp.deleteOnExit();
            outPath = temp.getAbsolutePath();

            InputStreamReader fReader = new InputStreamReader(GeneNetwork.class.getResourceAsStream("resources/post_stub.html"));
            reader = new BufferedReader(fReader);

            outputStream = new FileOutputStream(outPath);

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().equals("allthegraphmldatagoeshere")) {
                    writeEncodedString(Utilities.getString(this.createDocument()),
                            outputStream, true, true);
                } else {
                    outputStream.write((line + FileUtils.LINE_SEPARATOR).getBytes());
                    outputStream.flush();
                }
            }
            outputStream.flush();
        } catch (IOException e) {
            logger.error("Error writing cBio stub form to " + outPath);
            logger.error(e.getMessage());
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (outputStream != null) {
                outputStream.close();
            }
        }

        return outPath;
    }

    /**
     * Add the data specified by the score-types to our
     * network, using data from the tracks.
     * <p/>
     *
     * @param tracks
     * @param node_attributes
     */
    public void annotate(List<Track> tracks, Collection<String> node_attributes) {

        Set<Node> nodes = this.vertexSet();
        String name;
        ArrayList<RegionScoreType> types = new ArrayList<RegionScoreType>();
        for (String attr : node_attributes) {
            types.add(attribute_map.get(attr));
        }

        for (Node node : nodes) {
            name = getNodeKeyData(node, LABEL);

            ScoreData data = this.collectScoreData(name, tracks, types);


            //If we don't have any data to look at
            if (data != null) {
                //System.out.println("name: " + name + " total: " + data.getAvgScore() + " altered: " + data.getPercentAltered());
            } else {
                //System.out.println("name: " + name + " no data");
                continue;
            }

            float rel_data = data.getPercentAltered();
            if (rel_data == 0 && !Globals.isTesting()) {
                continue;
            }

            for (String attr : node_attributes) {
                Element newData = node.getOwnerDocument().createElement("data");
                newData.setAttribute(KEY, attr);
                newData.setTextContent("" + data.get(attribute_map.get(attr)));
                node.appendChild(newData);
            }

            //Set total
            Element newData = node.getOwnerDocument().createElement("data");
            newData.setAttribute(KEY, PERCENT_ALTERED);
            newData.setTextContent("" + data.getPercentAltered());

            node.appendChild(newData);
        }

        addSchema(Arrays.asList(PERCENT_ALTERED), "float", "node");
        addSchema(node_attributes, "float", "node");
    }

    public static ScoreData collectScoreData(String name, List<Track> tracks, Iterable<RegionScoreType> types) {
        int zoom = 0;

        List<NamedFeature> features = FeatureDB.getFeaturesList(name, Integer.MAX_VALUE);

        int numberSamples = features.size() * tracks.size();
        if (numberSamples == 0) {
            return null;
        }

        ScoreData<RegionScoreType, Float> results = new ScoreData(RegionScoreType.values().length);

        Set<String> anyAlteration = new HashSet<String>(numberSamples / 10);

        /*Track percentage of tracks have ANY of the specified alterations
        * Must account for the possibility of overlaps, can't double count
        */
        //int totalAnyAlteration = 0;

        for (RegionScoreType type : types) {

            //float totalScore = 0.0f;
            float percentAltered;
            int totalAltered = 0;

            for (NamedFeature feat : features) {
                int featStart = feat.getStart();
                int featEnd = feat.getEnd();
                for (Track track : tracks) {
                    //int anyChangeHere = 0;
                    float score = track.getRegionScore(feat.getChr(), featStart, featEnd, zoom,
                            type, Globals.isHeadless() ? null : FrameManager.getDefaultFrame().getName(), tracks);
                    //Note: Some methods return things like -Float.MAX_VALUE if they get confused
                    if (score > (-Float.MAX_VALUE + 1) && score > (Integer.MIN_VALUE + 1)) {
                        //totalScore += score;
                    } else {
                        continue;
                    }
                    totalAltered += score != 0.0 ? 1 : 0;
                    //anyChangeHere |= score != 0 ? 1 : 0;
                    if (score != 0.0) {
                        //Each track/feature pair represents a sample.
                        //All we require is a count of those altered
                        anyAlteration.add(feat.toString() + track.getName());
                    }
                }
                //totalAnyAlteration += anyChangeHere;
            }

            percentAltered = ((float) totalAltered) / numberSamples;
            //float avgScore = totalScore / numberSamples;
            results.put(type, percentAltered);
        }

        results.setPercentAltered(((float) anyAlteration.size()) / numberSamples);
        return results;

    }

    public void annotateAll(List<Track> tracks) {
        this.annotate(tracks, attribute_map.keySet());
    }


    public static class ScoreData<K, V> extends HashMap<K, V> {

        /**
         * Here we do not distinguish between any alteration value
         * or type. So 0,1,2,3 -> percentAltered = 3/4.
         * <p/>
         * Intended to represent the total fraction of samples with
         * ANY kind of alteration. So a sample that's mutated, amplified,
         * and upregulated would be counted once.
         */
        private float percentAltered;

        public ScoreData() {
        }

        public ScoreData(int size) {
            super(size);
        }


        public void setPercentAltered(float percentAltered) {
            this.percentAltered = percentAltered;
        }

        public float getPercentAltered() {
            return percentAltered;
        }

    }

    /**
     * Child class created so we can access the raw byte buffer,
     * and avoid making another copy.
     */
    private static class ByteArrayOutputStreamChild extends ByteArrayOutputStream {

        public ByteArrayOutputStreamChild() {
            super();
        }

        public ByteArrayOutputStreamChild(int size) {
            super(size);
        }

        public byte[] getBuf() {
            return this.buf;
        }
    }
}
