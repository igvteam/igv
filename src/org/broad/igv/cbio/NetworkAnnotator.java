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
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.tribble.Feature;
import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.*;
import java.net.URLEncoder;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for taking a network of genes and annotating them
 * with additional information.
 * <p/>
 * User: jacob
 * Date: 2012/02/09
 */
public class NetworkAnnotator {
    private static Logger logger = Logger.getLogger(NetworkAnnotator.class);

    private Document document;
    private Node graph;

    public static final String NODE_TAG = "node";
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


    static Map<String, RegionScoreType> attribute_map = new HashMap();

    static {
        attribute_map.put("PERCENT_MUTATED", RegionScoreType.MUTATION_COUNT);
        attribute_map.put("PERCENT_CNA_AMPLIFIED", RegionScoreType.AMPLIFICATION);
        attribute_map.put("PERCENT_CNA_HOMOZYGOUSLY_DELETED", RegionScoreType.DELETION);
    }

    public static ScoreData collectScoreData(String name, List<Track> tracks, RegionScoreType type) {
        int zoom;

        List<NamedFeature> features = FeatureDB.getFeaturesList(name, Integer.MAX_VALUE);

        int numberSamples = features.size() * tracks.size();
        if (numberSamples == 0) {
            return null;
        }
        float totalScore = 0.0f;
        float percentAltered;
        int totalAltered = 0;

        for (Feature feat : features) {
            for (Track track : tracks) {
                float score = track.getRegionScore(feat.getChr(), feat.getStart(), feat.getEnd(), zoom = -1,
                        type, Globals.isHeadless() ? null : FrameManager.getDefaultFrame().getName(), tracks);
                //Note: Some methods return things like -Float.MAX_VALUE if they get confused
                if (score > (-Float.MAX_VALUE + 1) && score > (Integer.MIN_VALUE + 1)) {
                    totalScore += score;
                }
                totalAltered += score != 0.0 ? 1 : 0;
            }
        }

        percentAltered = totalAltered / numberSamples;
        float avgScore = totalScore / numberSamples;
        return new ScoreData(avgScore, percentAltered);
    }

    public static NetworkAnnotator getFromCBIO(Iterable<String> geneList) {
        String query = HttpUtils.buildURLString(geneList, "+");
        String url = BASE_URL + "?" + GENE_LIST + "=" + query + "&" + common_parms;
        NetworkAnnotator annotator = new NetworkAnnotator();
        if (annotator.loadNetwork(url)) {
            return annotator;
        } else {
            return null;
        }
    }

    /**
     * Download and store gene network from cbio. Path can actually be to anywhere,
     * but it is expected to be graphml.
     *
     * @param path
     */
    public boolean loadNetwork(String path) {
        String error = null;
        try {
            InputStream cbioStream = ParsingUtils.openInputStreamGZ(new ResourceLocator(path));
            Document document = Utilities.createDOMDocumentFromXmlStream(cbioStream);
            this.document = document;
            this.graph = document.getElementsByTagName("graph").item(0);
        } catch (IOException e) {
            error = e.getMessage();
        } catch (ParserConfigurationException e) {
            error = e.getMessage();
        } catch (SAXException e) {
            error = e.getMessage();
        } finally {
            if (error != null) {
                logger.error(error);
                return false;
            } else {
                return true;
            }
        }
    }

    public NodeList getNodes() {
        return this.document.getElementsByTagName(NODE_TAG);
    }

    Document getDocument() {
        return this.document;
    }

    public void annotateAll(List<Track> tracks) {
        this.annotate(tracks, attribute_map.keySet());
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
            key = document.createElement("key");
            key.setAttribute("id", dataKey);
            //TODO id is supposed to unique, attr.name human readable.
            //Not quite sure of any reason they can't be the same.
            key.setAttribute("attr.name", dataKey);
            key.setAttribute("attr.type", dataType.toLowerCase());

            if (typeFor != null) {
                key.setAttribute("for", typeFor);
            }
            document.getFirstChild().appendChild(key);
        }
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

        NodeList nodes = getNodes();
        String name;
        Node node;
        for (int nn = 0; nn < nodes.getLength(); nn++) {
            node = nodes.item(nn);
            name = getNodeKeyData(node, LABEL);
            for (String attr : node_attributes) {
                RegionScoreType type = attribute_map.get(attr);
                ScoreData data = this.collectScoreData(name, tracks, type);
                //If we don't have any data to look at
                if (data == null) {
                    continue;
                }
                float rel_data = data.getPercentAltered();
                if (rel_data == 0 && !Globals.isTesting()) {
                    continue;
                }

                Element newData = this.document.createElement("data");
                newData.setAttribute(KEY, attr);
                newData.setTextContent("" + data.getPercentAltered());

                node.appendChild(newData);
            }
        }

        addSchema(node_attributes, "float", "node");
    }

    /**
     * The the value of a child node by the key.
     * Search is not recursive.
     * So getNodeKeyData(node, "key", "label") returns "JUN".
     * Note: This comment contains XML, it won't show
     * up properly when rendered.
     * <p/>
     * "<node id="3725"/>
     * <data key="label">JUN</data>
     * <data key="type">Protein</data>
     * <data key="RELATIONSHIP_XREF">HGNC:JUN;Entrez Gene:3725</data>
     * <data key="IN_QUERY">false</data>
     * </node>"
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
     * {@see getNodeAttrValue}
     *
     * @param node Node to search
     * @param key  Key to search for
     * @return String value of key found. null if not found
     */
    public static String getNodeKeyData(Node node, String key) {
        return getNodeAttrValue(node, KEY, key);
    }

    public int filterNodes(Predicate predicate) {
        NodeList nodeList = this.getNodes();
        int count = 0;
        for (int nn = 0; nn < nodeList.getLength(); nn++) {
            Node node = nodeList.item(nn);
            if (!predicate.evaluate(node)) {
                this.graph.removeChild(node);
                count++;
            }
        }
        return count;
    }

    boolean writeDocument(String outputFile) throws IOException {
        return writeDocument(document, outputFile);
    }

    /**
     * Write document to XML at outputFile. File is deleted if there
     * is an error writing out.
     *
     * @param document
     * @param outputFile
     * @return success
     * @throws java.io.IOException
     */
    public static boolean writeDocument(Document document, String outputFile) throws IOException {
        StreamResult streamResult;
        boolean success = false;
        try {

            TransformerFactory factory = TransformerFactory.newInstance();
            Transformer transformer = factory.newTransformer();

            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

            streamResult = new StreamResult(new StringWriter());
            DOMSource source = new DOMSource(document);
            transformer.transform(source, streamResult);
        } catch (TransformerException e) {
            return false;
        }

        String xmlString = streamResult.getWriter().toString();
        //System.out.println(xmlString);

        FileWriter fileWriter = new FileWriter(outputFile);
        try {
            fileWriter.write(xmlString);
            success = true;
        } catch (IOException e) {
            File outfile = new File(outputFile);
            if (outfile.exists()) {
                outfile.delete();
            }
        } finally {
            if (fileWriter != null) {
                fileWriter.close();
            }
        }
        return success;
    }


    public static class ScoreData {

        /**
         * The average of each score of a large number of tracks
         * So for 0,1,2,3 -> avgScore = 6/4
         */
        private float avgScore;

        /**
         * Here we do not distinguish between any alteration value.
         * So 0,1,2,3 -> percentAltered = 3/4.
         */
        private float percentAltered;

        public ScoreData(float avgScore, float percentAltered) {
            this.avgScore = avgScore;
            this.percentAltered = percentAltered;
        }


        public float getAvgScore() {
            return avgScore;
        }

        public float getPercentAltered() {
            return percentAltered;
        }

    }
}
