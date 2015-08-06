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

package org.broad.igv.das;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.tribble.CachingFeatureReader;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.w3c.dom.*;
import org.w3c.dom.traversal.DocumentTraversal;
import org.w3c.dom.traversal.NodeFilter;
import org.w3c.dom.traversal.TreeWalker;
import org.xml.sax.EntityResolver;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

//http://das.sanger.ac.uk/das/cosmic_mutations/features
//http://genome.ucsc.edu/cgi-bin/das/hg18/features?segment=1:1,100000;type=refGene

public class DASFeatureSource implements FeatureSource {

    private static Logger log = Logger.getLogger(DASFeatureSource.class);
    private static final WrappedIterator EMPTY__ITERATOR = new WrappedIterator(new ArrayList<Feature>().iterator());

    private String path;
    private String serverURL;
    private boolean isValid = true;
    private CachingFeatureReader reader;
    private int featureWindowSize = 1000000;   // 1 mbase window for loading features.  TODO -- need some dynamic way to set this
    private String type;
    private String[] parameters;
    Map<String, BasicFeature> groupFeatureCache = new HashMap(10000);


    public DASFeatureSource(ResourceLocator locator) throws MalformedURLException {
        URL url = new URL(locator.getPath());
        String host = url.getHost();
        String protocol = url.getProtocol();
        path = url.getPath();
        serverURL = protocol + "://" + host + path;

        String paramString = url.getQuery();
        if (paramString != null) {
            parameters = paramString.split(";");

            for (String param : parameters) {
                if (param.startsWith("type=")) {
                    type = param.substring(5);
                }
            }
        }

        if (locator.getPath().contains("genome.ucsc.edu") && type == null) {
            throw new DataLoadException("<html>Feature type is required for UCSC DAS tracks. <br>" +
                    "See http://www.broadinstitute.org/igv/LoadData for more details.",
                    locator.getPath());
        }

        //reader = new DasReader();
        reader = new CachingFeatureReader(new DasReader());
        reader.setBinSize(featureWindowSize);
    }

    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {
        return reader.query(chr, start, end);
    }

    public List<LocusScore> getCoverageScores(String chr, int i, int i1, int zoom) {
        return null;
    }

    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    public void setFeatureWindowSize(int size) {
        this.featureWindowSize = size;
        reader.setBinSize(size);
    }

    public Class getFeatureClass() {
        return BasicFeature.class;
    }


    public String getPath() {
        return path;
    }

    public String getType() {
        return type;
    }


    class DasReader implements FeatureReader {


        public CloseableTribbleIterator query(String chr, int start, int end, boolean contained) throws IOException {
            return query(chr, start, end);
        }

        public CloseableTribbleIterator query(String chr, int start, int end) throws IOException {

            int dasStart = start + 1;
            int dasEnd = end;

            if (isValid && !chr.equals("All")) {
                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                try {
                    groupFeatureCache.clear();
                    List<Feature> features = readFeatures(chr, dasStart, dasEnd);
                    if (features.size() < 1) {
                        return EMPTY__ITERATOR;
                    }
                    return new WrappedIterator(features.iterator());
                } finally {
                    groupFeatureCache.clear();
                    WaitCursorManager.removeWaitCursor(token);
                }
            }
            return EMPTY__ITERATOR;
        }

        // TODO Iterating not permitted -- throw exception?

        public CloseableTribbleIterator iterator() throws IOException {
            return null;
        }


        /**
         * Query for features in the given interval.  The coordinates are translated to DAS conventions.
         * <p/>
         * An end value <= 0 will result in a query for features over the whole chromosome.
         *
         * @param chr
         * @param start
         * @param end
         * @return
         */
        private List<Feature> readFeatures(String chr, int start, int end) {

            List<Feature> features = new ArrayList<Feature>();
            try {
                String dasChr = chr.startsWith("chr") ? chr.substring(3, chr.length()) : chr;

                String urlString = serverURL + "?" + "segment=" + dasChr;
                if (end > 0) urlString += ":" + start + "," + end;

                // Add parameters
                if (parameters != null) {
                    for (String param : parameters) {
                        if (!param.startsWith("segment")) {
                            urlString += ";" + param;
                        }
                    }
                }

                URL dataQuery = new URL(urlString);
                Document dom = getDocument(dataQuery);
                if (dom == null) {
                    return Collections.emptyList();
                }
                parseDocument(dom, chr, features);
                FeatureUtils.sortFeatureList(features);
                return features;

            } catch (IOException ioe) {
                throw new DataLoadException("Failed to reconnect with server", serverURL);
            }
        }


        private Document getDocument(URL query) {
            InputStream is = null;
            try {
                is = HttpUtils.getInstance().openConnectionStream(query);
                return createDocument(is);
            } catch (Exception e) {
                isValid = false;
                log.error(e);
                MessageUtils.showMessage("<html>The DAS Server: " + serverURL + " has returned an error or invalid data." +
                        "<br>" + e.getMessage());
                return null;
            } finally {
                if (is != null) {
                    try {
                        is.close();
                    } catch (IOException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }
        }


        private Document createDocument(InputStream inputStream)
                throws ParserConfigurationException, IOException, SAXException {
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            dbf.setValidating(false);
            DocumentBuilder documentBuilder = dbf.newDocumentBuilder();

            documentBuilder.setEntityResolver(new EntityResolver() {
                public InputSource resolveEntity(String publicId, String systemId) throws SAXException, IOException {
                    log.debug("Ignoring " + publicId + ", " + systemId);
                    return new InputSource(new StringReader(""));
                }
            });

            /*BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
            String nextLine;
            while((nextLine = br.readLine()) != null) {
                System.out.println("** " + nextLine);
            }
            */

            Document dom = documentBuilder.parse(inputStream);
            return dom;
        }

        public String getPath() {
            return serverURL;
        }

        private void parseDocument(Document dasDoc, String chr, List<Feature> features) {

            try {
                DocumentTraversal traversal = (DocumentTraversal) dasDoc;
                TreeWalker treewalker = traversal.createTreeWalker(
                        dasDoc.getDocumentElement(), NodeFilter.SHOW_ELEMENT, null, true);
                parseTree(treewalker, "FEATURE", chr, features);

            } catch (Exception ex) {
                log.error(ex);
                throw new DataLoadException("Error loading DAS resource (" + ex.toString() + ")", getPath());
            }
        }

        private List<Feature> parseTree(TreeWalker walker,
                                        String tag,
                                        String chr,
                                        List<Feature> features) {

            Node parent = walker.getCurrentNode();
            Element n = (Element) walker.firstChild();
            while (n != null) {
                if (n.getTagName().equalsIgnoreCase(tag)) {
                    Feature f = getFeature(walker, chr);
                    if (f != null) {
                        features.add(f);
                    }

                    n = (Element) walker.nextSibling();
                    continue;
                }
                parseTree(walker, tag, chr, features);
                n = (Element) walker.nextSibling();
            }
            walker.setCurrentNode(parent);
            return features;
        }

        private BasicFeature getFeature(TreeWalker walker, String chr) {

            String id;
            String label;

            String type = "";
            int start = 0;
            int end = 0;
            float score;
            Strand strand = Strand.NONE;

            int phase;
            String description = null;
            String link = null;
            String group = null;
            String groupLabel = null;
            String groupLink = null;


            Node featureNode = walker.getCurrentNode();

            //GET THE FEATURE ATTRIBUTES
            NamedNodeMap nnm = featureNode.getAttributes();
            Node tmpNode = nnm.getNamedItem("id");
            id = tmpNode.getTextContent();
            tmpNode = nnm.getNamedItem("label");
            label = tmpNode == null ? "" : tmpNode.getTextContent();

            for (Node n = walker.firstChild(); n != null;
                 n = walker.nextSibling()) {

                if (((Element) n).getTagName().equalsIgnoreCase("TYPE")) {
                    type = n.getTextContent();
                } else if (((Element) n).getTagName().equalsIgnoreCase("START")) {
                    start = Integer.parseInt(n.getTextContent()) - 1;
                } else if (((Element) n).getTagName().equalsIgnoreCase("END")) {
                    end = Math.max(start + 1, Integer.parseInt(n.getTextContent()));
                } else if (((Element) n).getTagName().equalsIgnoreCase("SCORE")) {
                    String scoreString = n.getTextContent();
                    if (!scoreString.equals("-")) score = Float.parseFloat(scoreString);
                } else if (((Element) n).getTagName().equalsIgnoreCase("PHASE")) {
                    String phaseString = n.getTextContent();
                    if (!phaseString.equals("-")) phase = Integer.parseInt(n.getTextContent());
                } else if (((Element) n).getTagName().equalsIgnoreCase("ORIENTATION")) {
                    String orientation = n.getTextContent();
                    if (orientation.equals("-")) {
                        strand = Strand.NEGATIVE;
                    } else if (orientation.equalsIgnoreCase("+")) {
                        strand = Strand.POSITIVE;
                    }
                } else if (((Element) n).getTagName().equalsIgnoreCase("NOTE")) {
                    if (description == null) {
                        description = "<html>" + n.getTextContent();
                    } else {
                        description += ("<br>" + n.getTextContent());
                    }
                } else if (((Element) n).getTagName().equalsIgnoreCase("GROUP")) {
                    nnm = n.getAttributes();
                    tmpNode = nnm.getNamedItem("id");
                    group = tmpNode.getTextContent();

                    tmpNode = nnm.getNamedItem("label");
                    if (tmpNode != null) {
                        groupLabel = tmpNode.getTextContent();
                    }

                    NodeList linkNodes = ((Element) n).getElementsByTagName("LINK");
                    if (linkNodes.getLength() > 0) {
                        Node ln = linkNodes.item(0);
                        Node hrefNode = ln.getAttributes().getNamedItem("href");
                        if (hrefNode != null) {
                            groupLink = hrefNode.getTextContent();
                        }
                    }
                    // TODO  group Note elements

                } else if (((Element) n).getTagName().equalsIgnoreCase("LINK")) {
                    NamedNodeMap tmpnnm = n.getAttributes();
                    Node tmpnode = tmpnnm.getNamedItem("href");
                    link = tmpnode.getTextContent();
                }

            }

            // TODO Rewind?  Why are we doing this?
            walker.setCurrentNode(featureNode);


            BasicFeature feature = null;
            if (group != null) {
                Exon exon = new Exon(chr, start, end, strand);

                feature = groupFeatureCache.get(group);
                if (feature == null) {
                    feature = new BasicFeature(exon.getChr(), exon.getStart(), exon.getEnd(), exon.getStrand());
                    feature.addExon(exon);
                    if (groupLink != null) {
                        feature.setURL(groupLink);
                    }
                    if (groupLabel != null) {
                        feature.setName(groupLabel);
                    } else {
                        feature.setName(label);
                    }
                    groupFeatureCache.put(group, feature);
                } else { // Seen before, just add the exon and exit
                    feature.addExon(exon);
                    return null;
                }

            } else {
                feature = new BasicFeature(chr, start, end);
                if (link != null) {
                    feature.setURL(link);
                }
                feature.setIdentifier(id);
                label = label.replace("?", " ");
                feature.setName(label);
                feature.setType(type);
                feature.setStrand(strand);
            }


            if (description != null) {
                description = description.replace("&nbsp;", "<br>");
                description = description.replace("<br><br><br>", "<br>");
                description = description.replace("<br><br>", "<br>");
                description = description.replace(":", ":&nbsp;");
                description = description.replace("?", " ");
                description = description + "<br>TYPE:&nbsp;" + type;
                feature.setDescription(description);
            }


            return feature;
        }

        public void close() throws IOException {
        }

        public List<String> getSequenceNames() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public Object getHeader() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

}
