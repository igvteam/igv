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
import org.apache.commons.collections.Predicate;
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

    private Map<String, Node> nodeTable;

    private List<Node> schema = new ArrayList<Node>();
    private NamedNodeMap graphAttr;
    private Document origDocument;


    static Map<String, RegionScoreType> attributeMap = new LinkedHashMap();

    /**
     * For each score type, we are only interested in the fraction of samples which are
     * modified in some way. We set thresholds here, in the form of [min max]
     */
    private static Map<String, float[]> bounds = new HashMap<String, float[]>(8);

    private Set<Node> rejectedNodes = new HashSet<Node>();
    private Set<Node> rejectedEdges = new HashSet<Node>();
    private boolean filtersFinalized;

    static {
        attributeMap.put("PERCENT_MUTATED", RegionScoreType.MUTATION_COUNT);
        attributeMap.put("PERCENT_CNA_AMPLIFIED", RegionScoreType.AMPLIFICATION);
        attributeMap.put("PERCENT_CNA_HOMOZYGOUSLY_DELETED", RegionScoreType.DELETION);
        attributeMap.put("PERCENT_MRNA_WAY_UP", RegionScoreType.EXPRESSION);
        attributeMap.put("PERCENT_MRNA_WAY_DOWN", RegionScoreType.EXPRESSION);
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
        bounds.put("PERCENT_MUTATED", new float[]{mut, max_val});

        //See GISTIC supplement, page 20
        float amp = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_AMPLIFICATION_THRESHOLD));
        bounds.put("PERCENT_CNA_AMPLIFIED", new float[]{amp, max_val});

        float del = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_DELETION_THRESHOLD));
        bounds.put("PERCENT_CNA_HOMOZYGOUSLY_DELETED", new float[]{del, max_val});

        //See GISTIC supplement, page 5, just gives greater than or less than 0
        float expUp = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_EXPRESSION_UP_THRESHOLD));
        bounds.put("PERCENT_MRNA_WAY_UP", new float[]{expUp, max_val});

        float expDown = Float.parseFloat(PreferenceManager.getInstance().get(PreferenceManager.CBIO_EXPRESSION_DOWN_THRESHOLD));
        bounds.put("PERCENT_MRNA_WAY_DOWN", new float[]{-max_val, -expDown});


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

    /**
     * Applies {@code predicate} to every element in {@code object}, and adds
     * any which return false to {@code rejectSet}. Intended for soft filtering.
     * There is an override here, where we never filter out a Node which
     * is marked as being part of the query
     *
     * @param predicate
     * @param objects
     * @param rejectSet
     * @return
     * @throws IllegalStateException If the filters have been finalized.
     */
    private int filter(Predicate predicate, Iterable<Node> objects, Set rejectSet) {
        if (filtersFinalized) {
            throw new IllegalStateException("Cannot filter after filtering has been finalized");
        }
        int filtered = 0;
        for (Node v : objects) {
            String in_query = getNodeAttrValue(v, KEY, "IN_QUERY");
            if (Boolean.parseBoolean(in_query)) {
                continue;
            }
            if (!predicate.evaluate(v)) {
                filtered += rejectSet.add(v) ? 1 : 0;
            }
        }
        return filtered;
    }

    public int filterNodes(Predicate predicate) {
        return this.filter(predicate, this.vertexSet(), rejectedNodes);
    }

    public int filterNodesRange(final String key, final float min, final float max) {
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

    public int filterEdges(Predicate predicate) {
        return this.filter(predicate, this.edgeSet(), rejectedEdges);
    }

    /**
     * Returns a set of nodes which have not been rejected by the filter
     *
     * @return
     */
    public Set<Node> vertexSetFiltered() {
        Set<Node> filteredSet = new HashSet<Node>(this.vertexSet());
        filteredSet.removeAll(rejectedNodes);
        return filteredSet;
    }

    /**
     * Return the set of edges connected to this node, whose neighbors
     * are not filtered out. If this node is filtered out, will
     * return the empty set.
     *
     * @param n
     * @return
     */
    public Set<Node> edgesOfFiltered(Node n) {
        if (rejectedNodes.contains(n)) {
            return new HashSet<Node>(0);
        }
        Set<Node> filteredEdges = new HashSet<Node>(this.edgesOf(n).size());
        for (Node edge : this.outgoingEdgesOf(n)) {
            if (!rejectedNodes.contains(this.getEdgeTarget(edge)) && !rejectedEdges.contains(edge)) {
                filteredEdges.add(edge);
            }
        }

        for (Node edge : this.incomingEdgesOf(n)) {
            if (!rejectedNodes.contains(this.getEdgeSource(edge)) && !rejectedEdges.contains(edge)) {
                filteredEdges.add(edge);
            }
        }
        return filteredEdges;
    }

    /**
     * Actually apply filters irrevocably, removing edges
     * and nodes from the graph
     *
     * @return True iff anything was removed due to filtering.
     *         Filters only get finalized if this call has an effect. So if
     *         this returns true, it means at least 1 node or edge has been removed,
     *         and the filters are finalized. If it returns false, nothing was removed
     *         and the graph is not finalized.
     */
    boolean finalizeFilters() {
        if (filtersFinalized) {
            throw new IllegalStateException("Cannot finalize filters twice");
        }
        boolean rejected = this.removeAllVertices(rejectedNodes);
        rejected |= this.removeAllEdges(rejectedEdges);
        this.clearAllFilters();
        this.filtersFinalized = rejected;
        return rejected;
    }


    public void clearNodeFilters() {
        if (filtersFinalized) {
            throw new IllegalStateException("Cannot clear filters after they have been finalized");
        }
        this.rejectedNodes.clear();
    }

    public void clearEdgeFilters() {
        if (filtersFinalized) {
            throw new IllegalStateException("Cannot clear filters after they have been finalized");
        }
        this.rejectedEdges.clear();
    }

    public void clearAllFilters() {
        this.clearNodeFilters();
        this.clearEdgeFilters();
    }

    public boolean pruneGraph() {
        Predicate min_connections = new Predicate<Node>() {
            public boolean evaluate(Node object) {
                return edgesOf(object).size() >= 1;
            }
        };
        return this.filterNodes(min_connections) > 0;
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
        String error = null;
        int numNodes = -1;
        try {
            InputStream cbioStream = ParsingUtils.openInputStreamGZ(new ResourceLocator(path));
            this.sourcePath = path;
            Document document = Utilities.createDOMDocumentFromXmlStream(cbioStream);

            //Cache the file
            if (HttpUtils.isURL(path)) {
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
            nodeTable = new HashMap<String, Node>(docNodes);
            for (int nn = 0; nn < docNodes; nn++) {
                Node node = nodes.item(nn);
                String label = node.getAttributes().getNamedItem("id").getTextContent();
                nodeTable.put(label, node);
                this.addVertex(node);
            }

            NodeList edges = document.getElementsByTagName(EDGE_TAG);
            int docEdges = edges.getLength();
            for (int ee = 0; ee < docEdges; ee++) {
                Node edge = edges.item(ee);
                NamedNodeMap attrs = edge.getAttributes();
                String source = attrs.getNamedItem("source").getTextContent();
                String target = attrs.getNamedItem("target").getTextContent();
                this.addEdge(nodeTable.get(source), nodeTable.get(target), edge);
            }

            numNodes = this.vertexSet().size();
        } catch (ParserConfigurationException e) {
            throw new IOException(e.getMessage());
        } catch (SAXException e) {
            throw new IOException(e.getMessage());
        }
        return numNodes;
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

    /**
     * Get the GraphML from this document. Filters will be finalized when this is called,
     * if they aren't already.
     *
     * @return
     */
    public Document createDocument() {
        if (!filtersFinalized) {
            this.finalizeFilters();
        }
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

    public void annotateAll(List<Track> tracks) {
        this.annotate(tracks, attributeMap.keySet());
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
