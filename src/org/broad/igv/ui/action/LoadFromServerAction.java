/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.ResourceTree;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXParseException;

import javax.swing.*;
import javax.xml.parsers.DocumentBuilderFactory;
import java.awt.event.ActionEvent;
import java.io.*;
import java.net.*;
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
        setToolTipText(UIConstants.LOAD_SERVER_DATA_TOOLTIP);
    }


    @Override
    public void actionPerformed(ActionEvent evt) {

        mainFrame.setStatusBarMessage("Loading ...");

        String urlString = PreferenceManager.getInstance().getDataServerURL();
        String genomeId = IGV.getInstance().getGenomeManager().getGenomeId();
        String genomeURL = urlString.replaceAll("\\$\\$", genomeId);
        try {

            // Add genome parameter for users with hardcoded reference to registry servlet
            if (IGVHttpUtils.isURL(genomeURL)) {
                genomeURL += "?genome=" + genomeId;
            }

            InputStream is = null;
            HttpURLConnection conn = null;
            LinkedHashSet<URL> urls = null;
            try {
                URL masterResourceFileURL = new URL(genomeURL);
                if (genomeURL.startsWith("ftp:")) {
                    is = IGVHttpUtils.openConnectionStream(masterResourceFileURL);
                } else if (IGVHttpUtils.isURL(genomeURL)) {
                    conn = (HttpURLConnection) IGVHttpUtils.openConnection(masterResourceFileURL);
                    conn.setReadTimeout(20000);
                    conn.setConnectTimeout(20000);
                    is = IGVHttpUtils.openHttpStream(masterResourceFileURL, conn);
                } else {
                    File file = new File(genomeURL.startsWith("file:") ? masterResourceFileURL.getFile() : genomeURL);
                    if (!file.exists()) {
                        MessageUtils.showMessage("ERROR: File does not exist: " + file.getAbsolutePath());
                        return;
                    }
                    is = new FileInputStream(file);
                }

                BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(is));
                urls = getResourceUrls(bufferedReader);
            } catch (IOException e) {
                if (conn.getResponseCode() == 401) {
                    String msg = "Authorization failure accessing resource: " + genomeURL;
                    MessageUtils.showMessage(msg);
                    log.error("Error accessing URL: " + genomeURL, e);
                }
            }
            finally {
                if (is != null) is.close();
                IGVHttpUtils.disconnect(conn);
            }

            if (urls == null || urls.isEmpty()) {
                JOptionPane.showMessageDialog(mainFrame.getMainFrame(),
                        "No datasets are available for the current genome (" + genomeId + ").");
            } else {

                List<ResourceLocator> locators = selectResources(urls);
                if (locators != null) {
                    mainFrame.loadTracks(locators);
                }
            }


        } catch (UnknownHostException e) {
            String msg = "The application could not contact the server: Unknown Host Failure\nHost: " + e.getMessage();
            log.error(msg);
            MessageUtils.showMessage(msg);

        } catch (ConnectException e) {
            String msg = "The application could not contact the server: Connection Failure!\n" + e.getMessage();
            log.error(msg);
            MessageUtils.showMessage(msg);
        }

        catch (Exception e) {
            String msg = "Could not load information from server!\n" + e.getMessage();
            log.error(msg);
            MessageUtils.showMessage(msg);

        } finally {
            mainFrame.showLoadedTrackCount();
        }

    }

    public List<ResourceLocator> selectResources(final LinkedHashSet<URL> xmlUrls) {

        if ((xmlUrls == null) || xmlUrls.isEmpty()) {
            log.error("No datasets are available from this server for the current genome (");
            return null;
        }

        StringBuffer buffer = new StringBuffer();
        boolean xmlParsingError = false;
        try {

            buffer.append("<html>The following urls could not be processed due to load failures:<br>");

            Document masterDocument = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();

            Element rootNode = masterDocument.createElement("Global");
            rootNode.setAttribute("name", "Available Datasets");
            rootNode.setAttribute("version", "1");

            masterDocument.appendChild(rootNode);

            // Merge all documents into one xml document for processing
            for (URL url : xmlUrls) {

                // Skip urls that have previously failed due to authorization
                if (failedURLs.contains(url.toString())) {
                    continue;
                }

                try {
                    HttpURLConnection connection = null;
                    InputStream is = null;
                    Document xmlDocument = null;
                    try {
                        connection = (HttpURLConnection) IGVHttpUtils.openConnection(url);
                        connection.setConnectTimeout(20000);
                        connection.setReadTimeout(20000);
                        is = IGVHttpUtils.openHttpStream(url, connection);
                        xmlDocument = Utilities.createDOMDocumentFromXmlStream(is);
                    }
                    catch (java.net.SocketTimeoutException e) {
                        xmlParsingError = true;
                        buffer.append("Error. Connection time out reading: " + url.toString());
                        continue;
                    }
                    catch (SAXParseException e) {
                        log.error("Invalid XML resource: " + url, e);

                        xmlParsingError = true;
                        buffer.append(url);
                        buffer.append("<br><i>");
                        if (url.toString().contains("iwww.broad")) {
                            buffer.append("File could not be loaded from the Broad Intranet");
                        } else {
                            buffer.append(e.getMessage());
                        }
                        buffer.append("");
                        continue;
                    } catch (FileNotFoundException e) {

                        String message = "Could not find file represented by " + url.toString();
                        log.error(message, e);

                        xmlParsingError = true;
                        buffer.append(url);
                        buffer.append("\t  [");
                        buffer.append(e.getMessage());
                        buffer.append("]\n");
                        continue;
                    }
                    catch (IOException e) {
                        String msg = "";
                        if (connection instanceof HttpURLConnection && ((HttpURLConnection) connection).getResponseCode() == 401) {
                            failedURLs.add(url.toString());
                            msg = "Authorization failure accessing resource: " + url;
                        } else {
                            msg = "Error accessing dataset list: " + e.toString();
                        }
                        MessageUtils.showMessage(msg);
                        log.error("Error accessing URL: " + url, e);
                    }
                    finally {
                        if (is != null) is.close();
                        IGVHttpUtils.disconnect(connection);
                    }

                    if (xmlDocument != null) {
                        NodeList elements = xmlDocument.getElementsByTagName("Global");
                        Element global = (Element) elements.item(0);
                        NodeList nodes = global.getChildNodes();
                        Element categoryNode = masterDocument.createElement("Category");
                        categoryNode.setAttribute("name", global.getAttribute("name"));
                        categoryNode.setAttribute("hyperlink", global.getAttribute("hyperlink"));
                        rootNode.appendChild(categoryNode);
                        int size = nodes.getLength();
                        for (int i = 0; i < size; i++) {
                            categoryNode.appendChild(masterDocument.importNode(nodes.item(i), true));
                        }
                    }
                } catch (Exception e) {
                    String message = "Cannot create an XML Document from " + url.toString();
                    log.error(message, e);
                    continue;
                }

            }
            if (xmlParsingError) {
                JOptionPane.showMessageDialog(mainFrame.getMainFrame(), buffer.toString());
            }


            /**
             * Resource Tree
             */
            HashSet<ResourceLocator> selectedLocators =
                    ResourceTree.getInstance().showResourceTreeDialog(mainFrame.getMainFrame(),
                            masterDocument, "Available Datasets");

            List<ResourceLocator> newLoadList = new ArrayList();

            Set<ResourceLocator> loadedResources = IGV.getInstance().getTrackManager().getDataResourceLocators();
            loadedResources.addAll(AttributeManager.getInstance().getLoadedResources());

            if (selectedLocators != null) {
                for (ResourceLocator locator : selectedLocators) {

                    // Don't reload data that is already loaded
                    if (loadedResources.contains(locator)) {
                        continue;
                    }

                    newLoadList.add(locator);
                }
            }

            return newLoadList;

        } catch (Exception e) {
            log.error("Could not load information from server", e);
            return null;
        } finally {
            if (xmlParsingError) {
                log.error(buffer.toString());
            }
        }
    }

    /**
     * Returns the complete list of URLs from the master registry file.
     *
     * @param bufferedReader
     * @return
     * @throws java.net.MalformedURLException
     * @throws java.io.IOException
     */
    private LinkedHashSet<URL> getResourceUrls
            (BufferedReader
                    bufferedReader)
            throws IOException {

        LinkedHashSet<URL> xmlFileUrls = new LinkedHashSet();
        while (true) {
            String xmlFileUrl = bufferedReader.readLine();
            if ((xmlFileUrl == null) || (xmlFileUrl.trim().length() == 0)) {
                break;
            }
            xmlFileUrl = xmlFileUrl.trim();
            xmlFileUrls.add(new URL(xmlFileUrl));
        }

        return xmlFileUrls;
    }


}
