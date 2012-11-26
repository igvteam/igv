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

package org.broad.igv.cbio;

import biz.source_code.base64Coder.Base64Coder;
import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.Collections2;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.*;
import org.jgrapht.EdgeFactory;
import org.jgrapht.graph.DirectedMultigraph;
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
 * <p/>
 * We implement soft filtering, where filter methods remove
 * any node or edge for which predicate returns false. All
 * original data can be restored with {@link #reset}, which resets
 * to the last time loadNetwork was called
 * User: jacob
 * Date: 2012/01/31
 */
public class GeneNetwork extends DirectedMultigraph<Node, Node> {

    private Logger log = Logger.getLogger(GeneNetwork.class);

    public static final String NODE_TAG = "node";
    public static final String EDGE_TAG = "edge";
    public static final String KEY = "key";
    public static final String LABEL = "label";


    /**
     * URL that cbio will use when service is released
     */
    static final String REAL_URL = "http://www.cbioportal.org/public-portal/network.do";
    /**
     * URL they use for testing
     */
    static final String TEST_URL = "http://awabi.cbio.mskcc.org/public-portal/network.do";
    static String BASE_URL = REAL_URL;
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

    private List<Node> schema = new ArrayList<Node>();
    private NamedNodeMap graphAttr;


    private Document origDocument;


    static Map<String, RegionScoreType> attributeMap = new LinkedHashMap();

    /**
     * For each score type, we are only interested in the fraction of samples which are
     * modified in some way. We set thresholds here, in the form of [min max]
     */
    private static Map<String, float[]> bounds = new HashMap<String, float[]>(8);

    /**
     * We keep the original graph so we can reset filters
     */
    private DirectedMultigraph<Node, Node> origGraph;

    /**
     * We only filter certain node types (e.g. "Protein", others we just leave alone)
     */
    private static final Set<String> geneTypes = new HashSet<String>(5);

    static final Predicate<Node> isGene;
    static final Predicate<Node> isNotGene;
    static final Predicate<Node> inQuery;

    static final String PERCENT_MUTATED = "PERCENT_MUTATED";
    static final String PERCENT_CNA_AMPLIFIED = "PERCENT_CNA_AMPLIFIED";
    static final String PERCENT_CNA_HOMOZYGOUSLY_DELETED = "PERCENT_CNA_HOMOZYGOUSLY_DELETED";
    static final String PERCENT_MRNA_WAY_UP = "PERCENT_MRNA_WAY_UP";
    static final String PERCENT_MRNA_WAY_DOWN = "PERCENT_MRNA_WAY_DOWN";

    static {
        attributeMap.put(PERCENT_MUTATED, RegionScoreType.MUTATION_COUNT);
        attributeMap.put(PERCENT_CNA_AMPLIFIED, RegionScoreType.AMPLIFICATION);
        attributeMap.put(PERCENT_CNA_HOMOZYGOUSLY_DELETED, RegionScoreType.DELETION);
        attributeMap.put(PERCENT_MRNA_WAY_UP, RegionScoreType.EXPRESSION);
        attributeMap.put(PERCENT_MRNA_WAY_DOWN, RegionScoreType.EXPRESSION);

        geneTypes.add("Protein");
        isGene = new Predicate<Node>() {
            @Override
            public boolean apply(Node object) {
                String type = getNodeKeyData(object, "TYPE");
                return geneTypes.contains(type);
            }
        };

        isNotGene = Predicates.not(isGene);

        inQuery = new Predicate<Node>() {

            @Override
            public boolean apply(Node object) {
                String in_query = getNodeAttrValue(object, KEY, "IN_QUERY");
                return Boolean.parseBoolean(in_query);
            }
        };
    }

    /**
     * The thresholds determine whether a particular sample is altered in that fashion.
     * The number is calculated using track.getRegionScore
     *
     * Mutation will be the number of mutations in a given gene (integer)
     *
     * For the others, we take the numerical value that has been read from the data source,
     * and average it over the region of each gene.
     *
     */
    static {
        float max_val = 2 << 10;

        float mut = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_MUTATION_THRESHOLD));
        bounds.put(PERCENT_MUTATED, new float[]{mut, max_val});

        //See GISTIC supplement, page 20
        float amp = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_AMPLIFICATION_THRESHOLD));
        bounds.put(PERCENT_CNA_AMPLIFIED, new float[]{amp, max_val});

        float del = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_DELETION_THRESHOLD));
        bounds.put(PERCENT_CNA_HOMOZYGOUSLY_DELETED, new float[]{del, max_val});

        //See GISTIC supplement, page 5, just gives greater than or less than 0
        float expUp = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_EXPRESSION_UP_THRESHOLD));
        bounds.put(PERCENT_MRNA_WAY_UP, new float[]{expUp, max_val});

        float expDown = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_EXPRESSION_DOWN_THRESHOLD));
        bounds.put(PERCENT_MRNA_WAY_DOWN, new float[]{-max_val, -expDown});


    }

    public static final String PERCENT_ALTERED = "PERCENT_ALTERED";


    private String sourcePath;

    /**
     * Generally for testing. The path from which data was loaded.
     *
     * @return
     */
    String getSourcePath() {
        return sourcePath;
    }

    GeneNetwork() {
        this(Node.class);
    }

    GeneNetwork(EdgeFactory edgeFactory) {
        super(edgeFactory);
    }

    private GeneNetwork(Class clazz) {
        super(clazz);
    }

    /**
     * Applies {@code predicate} to every element in {@code object}, and adds
     * any which return false to {@code rejectSet}. Intended for soft filtering.
     *
     * @param predicate
     * @param objects
     * @return
     */
    private Set<Node> filter(Predicate<Node> predicate, Collection<Node> objects) {
        Set<Node> rejectedSet = new HashSet<Node>(objects.size());
        for (Node v : objects) {
            if (!predicate.apply(v)) {
                rejectedSet.add(v);
            }
        }
        return rejectedSet;
    }

    private int filterNodes(Predicate<Node> predicate) {
        Set<Node> rejectedSet = this.filter(predicate, this.vertexSet());
        this.removeAllVertices(rejectedSet);
        return rejectedSet.size();
    }

    /**
     * Applies this predicate only to the genes. Any nodes
     * which are NOT genes are automatically kept.
     * <p/>
     * There is an override here, where we never filter out a Node which
     * is marked as being part of the query
     *
     * @param predicate
     * @return
     */
    public int filterGenes(Predicate<Node> predicate) {
        Predicate<Node> genePredicate = Predicates.or(predicate, isNotGene);
        Predicate<Node> finalPredicate = Predicates.or(genePredicate, inQuery);
        return this.filterNodes(finalPredicate);
    }

    public int filterGenesRange(final String key, final float min, final float max) {
        Predicate<Node> pred = new Predicate<Node>() {
            @Override
            public boolean apply(Node object) {
                String sval = getNodeKeyData(object, key);
                if (sval != null) {
                    float fval = Float.parseFloat(sval);
                    return fval >= min && fval <= max;
                }
                return false;
            }
        };

        return this.filterGenes(pred);
    }

    public int filterEdges(Predicate<Node> predicate) {
        Set<Node> rejectedSet = this.filter(predicate, this.edgeSet());
        this.removeAllEdges(rejectedSet);
        return rejectedSet.size();
    }

    /**
     * Returns a set of gene nodes which have not been rejected by the filter
     *
     * @return
     */
    public Collection<Node> geneVertexes() {
        return Collections2.filter(vertexSet(), isGene);
    }

    /**
     * Reset the graph to the Nodes which were originally
     * contained in the {@link #loadNetwork(String)} call.
     * Since we do shallow copying, if the nodes themselves have been
     * altered, those alterations will be preserved.
     */
    public void reset() {
        if (this.origGraph == null) {
            throw new IllegalStateException("Have no original graph to which to reset");
        }

        this.removeAllVertices(new HashSet<Node>(this.vertexSet()));
        this.removeAllEdges(new HashSet<Node>(this.edgeSet()));

        copyGraph(this.origGraph);
    }

    public boolean pruneGraph() {
        Predicate<Node> min_connections = new Predicate<Node>() {
            public boolean apply(Node object) {
                return edgesOf(object).size() >= 1;
            }
        };
        return this.filterGenes(min_connections) > 0;
    }

    /**
     * Hash the url to get a file name, locate it in temp directory
     *
     * @param url
     * @return
     */
    static File getCachedFile(String url) {
        String cachedFileName = Math.abs(url.hashCode()) + "_tmp.xml";
        File cachedFile = new File(DirectoryManager.getCacheDirectory(), cachedFileName);
        return cachedFile;
    }

    static String getURLForGeneList(Iterable<String> geneList) {
        String query = HttpUtils.buildURLString(geneList, "+");
        String url = BASE_URL + "?" + GENE_LIST + "=" + query + "&" + common_parms;
        return url;
    }

    public static GeneNetwork getFromCBIO(Iterable<String> geneList) throws IOException {
        String url = getURLForGeneList(geneList);

        File cachedFile = getCachedFile(url);
        if (cachedFile.exists()) {
            url = cachedFile.getAbsolutePath();
        }

        GeneNetwork network = new GeneNetwork();
        network.loadNetwork(url);
        return network;
    }

    /**
     * Load graphml data from provide path.
     *
     * @param path
     * @throws IOException If something went wrong loading data.
     *                     Note we wrap other exception types in this.
     */
    public int loadNetwork(String path) throws IOException {
        origGraph = new DirectedMultigraph<Node, Node>(Node.class);
        String error = null;
        int numNodes = -1;
        try {
            InputStream cbioStream = ParsingUtils.openInputStreamGZ(new ResourceLocator(path));
            this.sourcePath = path;
            Document document = Utilities.createDOMDocumentFromXmlStream(cbioStream);

            //Cache the file
            if (HttpUtils.isRemoteURL(path)) {
                File cacheFile = getCachedFile(path);
                try {
                    this.exportDocument(document, cacheFile.getAbsolutePath());
                } catch (IOException e) {
                    //No biggy, we just don't cache the file
                    log.error("Error caching file: " + e);
                    cacheFile.delete();
                }
                cacheFile.deleteOnExit();
            }

            this.origDocument = document;

            //Read schema from top and save it
            addToSchema(document.getElementsByTagName("key"));

            graphAttr = document.getElementsByTagName("graph").item(0).getAttributes();


            NodeList nodes = document.getElementsByTagName(NODE_TAG);
            //Generate the graph itself. First add the nodes, then the edges
            int docNodes = nodes.getLength();
            Map<String, Node> nodeTable = new HashMap<String, Node>(docNodes);
            for (int nn = 0; nn < docNodes; nn++) {
                Node node = nodes.item(nn);
                String label = node.getAttributes().getNamedItem("id").getTextContent();
                nodeTable.put(label, node);
                origGraph.addVertex(node);
            }

            NodeList edges = document.getElementsByTagName(EDGE_TAG);
            int docEdges = edges.getLength();
            for (int ee = 0; ee < docEdges; ee++) {
                Node edge = edges.item(ee);
                NamedNodeMap attrs = edge.getAttributes();
                String source = attrs.getNamedItem("source").getTextContent();
                String target = attrs.getNamedItem("target").getTextContent();
                origGraph.addEdge(nodeTable.get(source), nodeTable.get(target), edge);
            }

            this.reset();
            numNodes = this.vertexSet().size();
        } catch (ParserConfigurationException e) {
            throw new IOException(e);
        } catch (SAXException e) {
            throw new IOException(e);
        }
        return numNodes;
    }

    /**
     * Shallow-copy the graph from origGraph to this.
     * Does NOT remove anything from this graph, or alter origGraph
     *
     * @param origGraph
     */
    private void copyGraph(DirectedMultigraph<Node, Node> origGraph) {
        for (Node sourceV : origGraph.vertexSet()) {
            this.addVertex(sourceV);
            for (Node edge : origGraph.outgoingEdgesOf(sourceV)) {
                Node targetV = origGraph.getEdgeTarget(edge);
                if (!this.containsVertex(targetV)) {
                    this.addVertex(targetV);
                }
                this.addEdge(sourceV, targetV, edge);
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
                String nodeValue = label.getNodeValue();
                if (nodeValue.compareToIgnoreCase(attrValue) == 0) {
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

    /**
     * Get the GraphML from this document, keeping only Nodes which passed filtering.
     *
     * @return
     */
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
            throw new RuntimeException("Error outputting graph", e);
        }
    }

    public static int writeEncodedString(String string, OutputStream outputStream, boolean gzip, boolean base64encode) throws IOException {

        byte[] byteData;

        if (gzip) {
            ByteArrayOutputStream gmlByteStream = new ByteArrayOutputStream(string.length() / 20);
            GZIPOutputStream gzipOutputStream = new GZIPOutputStream(gmlByteStream);
            gzipOutputStream.write(string.getBytes());
            gzipOutputStream.finish();
            byteData = gmlByteStream.toByteArray();
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
     * @param document   The Document to write out
     * @param outputPath
     * @return success
     * @throws java.io.IOException
     */
    private int exportDocument(Document document, String outputPath) throws IOException {
        boolean gzip = outputPath.endsWith(".gz");

        String xmlString = Utilities.getString(document);

        OutputStream stream = new FileOutputStream(outputPath);
        int count = writeEncodedString(xmlString, stream, gzip, false);

        stream.flush();
        stream.close();

        return count;
    }

    /**
     * Calls {@link #exportDocument} with the Document created from {@link #createDocument}
     *
     * @param outputFile
     * @return
     * @throws IOException
     */
    int exportGraph(String outputFile) throws IOException {
        return exportDocument(this.createDocument(), outputFile);
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
            log.error("Error writing cBio stub form to " + outPath);
            log.error(e.getMessage());
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

    public void annotateAll(List<Track> tracks) {
        this.annotate(tracks, attributeMap.keySet());
    }

    /**
     * Add the data specified by the score-types to our
     * network, using data from the tracks.
     * <p/>
     *
     * @param tracks
     * @param nodeAttributes
     */
    public void annotate(List<Track> tracks, Collection<String> nodeAttributes) {
        Set<Node> nodes = this.vertexSet();
        String name;

        for (Node node : nodes) {
            name = getNodeKeyData(node, LABEL);

            ScoreData data = this.collectScoreData(name, tracks, nodeAttributes);

            //Don't add annotation if gene has no alteration?
            float relData = data.getPercentAltered();
            if (relData == 0 && !Globals.isTesting()) {
                continue;
            }

            for (String attr : nodeAttributes) {
                Element newData = node.getOwnerDocument().createElement("data");
                newData.setAttribute(KEY, attr);
                newData.setTextContent("" + data.get(attr));
                node.appendChild(newData);
            }

            //Set total
            Element newData = node.getOwnerDocument().createElement("data");
            newData.setAttribute(KEY, PERCENT_ALTERED);
            newData.setTextContent("" + data.getPercentAltered());

            node.appendChild(newData);
        }

        addSchema(Arrays.asList(PERCENT_ALTERED), "float", "node");
        addSchema(nodeAttributes, "float", "node");
    }

    public ScoreData collectScoreData(String name, List<Track> tracks, Iterable<String> attributes) {
        int zoom = 0;

        List<NamedFeature> features = FeatureDB.getFeaturesList(name, Integer.MAX_VALUE);

        //If we are viewing a gene list, use the frame
        List<ReferenceFrame> frames = Globals.isHeadless() ? null : FrameManager.getFrames();
        ReferenceFrame frame = Globals.isHeadless() ? null : FrameManager.getDefaultFrame();
        if (frames != null) {
            for (ReferenceFrame frm : frames) {
                if (frm.getName().equalsIgnoreCase(name)) {
                    frame = frm;
                }
            }
        }
        String frameName = frame != null ? frame.getName() : null;

        ScoreData<String, Float> results = new ScoreData(RegionScoreType.values().length);

        int initCapacity = tracks.size() / 10;

        //The names of all samples these tracks cover
        Set<String> allSamples = new HashSet<String>(initCapacity);

        //Each track/feature pair represents a region of a sample.
        //We store whether that sample has been altered in ANY way
        Set<String> anyAlteration = new HashSet<String>(initCapacity);


        //Set of samples which have data for this type
        Set<String> samplesForType = new HashSet<String>(initCapacity);
        //Set of samples which have been altered, using this type.
        Set<String> alteredSamplesForType = new HashSet<String>(initCapacity);

        for (String attr : attributes) {
            if (!bounds.containsKey(attr)) {
                throw new IllegalArgumentException("Have no bounds for " + attr);
            }

            RegionScoreType type = attributeMap.get(attr);
            float[] curBounds = bounds.get(attr);

            samplesForType.clear();
            alteredSamplesForType.clear();


            for (NamedFeature feat : features) {
                if (!name.equalsIgnoreCase(feat.getName())) {
                    continue;
                }
                int featStart = feat.getStart();
                int featEnd = feat.getEnd();
                for (Track track : tracks) {
                    if (!track.isVisible()) {
                        continue;
                    }
                    String sample = track.getSample();

                    //If track is wrong type, or if sample has already been marked altered,
                    //no further information can be gained
                    if (!track.isRegionScoreType(type) || alteredSamplesForType.contains(sample)) {
                        //if(alteredSamplesForType.contains(sample)) assert samplesForType.contains(sample);
                        continue;
                    }

                    samplesForType.add(sample);

                    float score = track.getRegionScore(feat.getChr(), featStart, featEnd, zoom,
                            type, frameName, tracks);

                    if (score >= curBounds[0] && score <= curBounds[1] && !Float.isNaN(score)) {
                        alteredSamplesForType.add(sample);
                    }
                }
            }

            allSamples.addAll(samplesForType);
            anyAlteration.addAll(alteredSamplesForType);

            float fractionAltered = ((float) alteredSamplesForType.size()) / samplesForType.size();
            results.put(attr, fractionAltered);
        }

        results.setPercentAltered(((float) anyAlteration.size()) / allSamples.size());
        return results;

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

}
