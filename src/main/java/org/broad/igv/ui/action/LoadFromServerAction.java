/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.logging.*;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.ResourceTree;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * @author jrobinso
 */
public class LoadFromServerAction extends MenuAction {

    static Logger log = LogManager.getLogger(LoadFromServerAction.class);

    IGV mainFrame;

    static Map<String, Set<String>> cachedNodeURLs = new HashMap();

    // Keep track of authorization failures so user isn't constantly harranged
    static HashSet<String> failedURLs = new HashSet();


    public LoadFromServerAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    public static String getGenomeDataURL(String genomeId) {
        String urlString = PreferencesManager.getPreferences().getDataServerURL();
        String genomeURL = urlString == null ? null : urlString.replaceAll("\\$\\$", genomeId);
        return genomeURL;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {

        mainFrame.setStatusBarMessage("Loading ...");
        String genomeId = GenomeManager.getInstance().getGenomeId();

        try {

            Set<String> nodeURLs = getNodeURLs(genomeId);

            if (nodeURLs == null || nodeURLs.isEmpty()) {
                MessageUtils.showMessage("No datasets are available for the current genome (" + genomeId + ").");
            } else {

                List<ResourceLocator> locators = loadNodes(nodeURLs);
                if (locators != null) {
                    mainFrame.loadTracks(locators);
                }
            }


        } finally {
            mainFrame.showLoadedTrackCount();
        }

    }

    public static Set<String> getNodeURLs(String genomeId) {

        String genomeURL = getGenomeDataURL(genomeId);

        Set<String> nodeURLs;
        synchronized (cachedNodeURLs) {
            nodeURLs = cachedNodeURLs.get(genomeId);
            if (nodeURLs == null) {
                if (genomeURL != null && genomeURL.length() > 0) {
                    InputStream is = null;
                    try {
                        is = ParsingUtils.openInputStreamGZ(new ResourceLocator(genomeURL));
                        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(is));
                        nodeURLs = getResourceUrls(bufferedReader);
                    } catch (IOException e) {
                        // This is common and not an error, a load-from-server "data registry" file is optional
                        nodeURLs = Collections.emptySet();
                    } finally {
                        if (is != null) {
                            try {
                                is.close();
                            } catch (IOException e) {
                                log.error("Error closing input stream", e);
                            }
                        }
                    }
                    cachedNodeURLs.put(genomeId, nodeURLs);
                }
            }
        }

        return nodeURLs;
    }

    private List<ResourceLocator> loadNodes(final Set<String> xmlUrls) {

        if ((xmlUrls == null) || xmlUrls.isEmpty()) {
            log.error("No datasets are available from this server for the current genome (");
            return null;
        }

        try {

            Document masterDocument = createMasterDocument(xmlUrls);

            /**
             * Resource Tree
             */
            LinkedHashSet<ResourceLocator> selectedLocators =
                    ResourceTree.showResourceTreeDialog(mainFrame.getMainFrame(),
                            masterDocument, "Available Datasets");

            List<ResourceLocator> newLoadList = new ArrayList();

            if (selectedLocators != null) {
                for (ResourceLocator locator : selectedLocators) {

                    // Don't reload data that is already loaded
                    if (IGV.getInstance().getDataResourceLocators().contains(locator)) {
                        continue;
                    }

                    newLoadList.add(locator);
                }
            }

            return newLoadList;

        } catch (Exception e) {
            log.error("Could not load information from server", e);
            return null;
        }
    }


    public static Document createMasterDocument(Collection<String> xmlUrls) throws ParserConfigurationException {

        StringBuffer buffer = new StringBuffer();

        Document masterDocument = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();

        Element rootNode = masterDocument.createElement("Global");
        rootNode.setAttribute("name", "Available Datasets");
        rootNode.setAttribute("version", "1");

        masterDocument.appendChild(rootNode);


        // Merge all documents into one xml document for processing
        for (String url : xmlUrls) {

            // Skip urls that have previously failed due to authorization
            if (failedURLs.contains(url)) {
                continue;
            }

            try {
                Document xmlDocument = readXMLDocument(url, buffer);

                if (xmlDocument != null) {
                    Element global = xmlDocument.getDocumentElement();
                    masterDocument.getDocumentElement().appendChild(masterDocument.importNode(global, true));

                }
            } catch (Exception e) {
                String message = "Cannot create an XML Document from " + url.toString();
                log.error(message, e);
                continue;
            }

        }
        if (buffer.length() > 0) {
            String message = "<html>The following urls could not be processed due to load failures:<br>" + buffer.toString();
            MessageUtils.showMessage(message);
        }

        return masterDocument;

    }

    private static Document readXMLDocument(String url, StringBuffer errors) {
        InputStream is = null;
        Document xmlDocument = null;
        //long t0 = System.currentTimeMillis();
        try {
            is = ParsingUtils.openInputStreamGZ(new ResourceLocator(url));
            xmlDocument = Utilities.createDOMDocumentFromXmlStream(is);
            xmlDocument = resolveIncludes(xmlDocument, errors);

        } catch (SAXException e) {
            log.error("Invalid XML resource: " + url, e);
            errors.append(url + "<br><i>" + e.getMessage());
        } catch (java.net.SocketTimeoutException e) {
            log.error("Connection time out", e);
            errors.append(url + "<br><i>Connection time out");
        } catch (IOException e) {
            log.error("Error accessing " + url.toString(), e);
            errors.append(url + "<br><i>" + e.getMessage());
        } catch (ParserConfigurationException e) {
            log.error("Parser configuration error for:" + url, e);
            errors.append(url + "<br><i>" + e.getMessage());
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    log.error("Error closing stream for: " + url, e);
                }
            }
        }
        //double dt = (System.currentTimeMillis() - t0) / 1000.0;
        //System.out.println(url + "  " + dt);
        return xmlDocument;
    }

    private static Document resolveIncludes(Document document, StringBuffer errors) {

        NodeList includeNodes = document.getElementsByTagName("Include");
        if (includeNodes.getLength() == 0) {
            return document;
        }

        int size = includeNodes.getLength();
        // Copy the nodes as we'll be modifying the tree.  This is neccessary!
        Node[] tmp = new Node[size];
        for (int i = 0; i < size; i++) {
            tmp[i] = includeNodes.item(i);
        }

        for (Node item : tmp) {
            NamedNodeMap nodeMap = item.getAttributes();
            if (nodeMap == null) {
                log.warn("XML node " + item.getNodeName() + " has no attributes");
            } else {
                Attr path = (Attr) item.getAttributes().getNamedItem("path");
                if (path == null) {
                    log.warn("XML node " + item.getNodeName() + " is missing a path attribute");
                } else {
                    Node parent = item.getParentNode();

                    //log.warn("Loading node " + path.getValue());
                    Document doc = readXMLDocument(path.getValue(), errors);
                    if (doc != null) {
                        Element global = doc.getDocumentElement();
                        Node expandedNode = parent.getOwnerDocument().importNode(global, true);
                        parent.replaceChild(expandedNode, item);
                    }
                }
            }
        }


        return document;

    }

    /**
     * Returns the complete list of URLs from the master registry file.
     *
     * @param bufferedReader
     * @return
     * @throws java.net.MalformedURLException
     * @throws java.io.IOException
     */
    private static LinkedHashSet<String> getResourceUrls(BufferedReader bufferedReader)
            throws IOException {

        LinkedHashSet<String> xmlFileUrls = new LinkedHashSet();
        String xmlFileUrl;
        while ((xmlFileUrl = bufferedReader.readLine()) != null) {

            xmlFileUrl = xmlFileUrl.trim();

            if (xmlFileUrl.length() == 0 || xmlFileUrl.startsWith("#")) {
                continue;
            }

            xmlFileUrls.add(xmlFileUrl);
        }

        return xmlFileUrls;
    }


}
