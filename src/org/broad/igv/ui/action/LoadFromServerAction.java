/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.batik.util.gui.xmleditor.XMLDocument;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.ResourceTree;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.*;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

import javax.swing.*;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 */
public class LoadFromServerAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFromServerAction.class);
    IGV mainFrame;

    // Keep track of authorization failures so user isn't constantly harranged
    static HashSet<String> failedURLs = new HashSet();


    public LoadFromServerAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    public static String getGenomeDataURL(String genomeId) {
        String urlString = PreferenceManager.getInstance().getDataServerURL();
        String genomeURL = urlString.replaceAll("\\$\\$", genomeId);
        return genomeURL;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {

        mainFrame.setStatusBarMessage("Loading ...");
        String genomeId = GenomeManager.getInstance().getGenomeId();
        String genomeURL = getGenomeDataURL(genomeId);

        try {

            LinkedHashSet<String> nodeURLs = getNodeURLs(genomeURL);

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

    public static LinkedHashSet<String> getNodeURLs(String genomeURL) {

        InputStream is = null;
        LinkedHashSet<String> nodeURLs = null;
        try {
            is = ParsingUtils.openInputStreamGZ(new ResourceLocator(genomeURL));
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(is));
            nodeURLs = getResourceUrls(bufferedReader);
        } catch (IOException e) {
            MessageUtils.showMessage("Error loading the data registry file: " + e.toString());
            log.error("Error loading genome registry file", e);
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    log.error("Error closing input stream", e);
                }
            }
        }

        return nodeURLs;
    }

    public List<ResourceLocator> loadNodes(final LinkedHashSet<String> xmlUrls) {

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
                    ResourceTree.getInstance().showResourceTreeDialog(mainFrame.getMainFrame(),
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


    private Document createMasterDocument(Collection<String> xmlUrls) throws ParserConfigurationException {

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

    private Document readXMLDocument(String url, StringBuffer errors) {
        InputStream is = null;
        Document xmlDocument = null;
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
        return xmlDocument;
    }

    private Document resolveIncludes(Document document, StringBuffer errors) {

        NodeList includeNodes = document.getElementsByTagName("Include");
        if (includeNodes.getLength() == 0) {
            return document;
        }

        int size = includeNodes.getLength();
        // Copy the nodes as we'll be modifying the tree.  This is neccessary!
        Node [] tmp = new Node[size];
        for(int i=0; i<size; i++) {
            tmp[i] = includeNodes.item(i);
        }

        for (Node item : tmp) {
            NamedNodeMap nodeMap = item.getAttributes();
            if (nodeMap == null) {
                log.info("XML node " + item.getNodeName() + " has no attributes");
            } else {
                Attr path = (Attr) item.getAttributes().getNamedItem("path");
                if (path == null) {
                    log.info("XML node " + item.getNodeName() + " is missing a path attribute");
                } else {
                    Node parent = item.getParentNode();

                    //log.info("Loading node " + path.getValue());
                    Document doc = this.readXMLDocument(path.getValue(), errors);
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
        while (true) {
            String xmlFileUrl = bufferedReader.readLine();
            if ((xmlFileUrl == null) || (xmlFileUrl.trim().length() == 0)) {
                break;
            }
            xmlFileUrl = xmlFileUrl.trim();
            xmlFileUrls.add(xmlFileUrl);
        }

        return xmlFileUrls;
    }


}
