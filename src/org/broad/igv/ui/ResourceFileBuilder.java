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
package org.broad.igv.ui;

//~--- non-JDK imports --------------------------------------------------------


import org.broad.igv.ui.color.ColorUtilities;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author eflakes
 */
public class ResourceFileBuilder {

    private static Document masterDocument = null;

    // private static Element rootCategoryNode = null;
    private static String DEFAULT_SERVER_URL = "http://www.broadinstitute.org/dataserver/data";

    private static String DEFAULT_OUTPUT_FILE = "tempResourceFile.xml";

    boolean epigenetics;

    Element rootNode = null;

    /**
     * Method description
     *
     * @param inputDir
     */
    public void process(File inputDir) {
        try {
            (new ResourceFileBuilder()).process(inputDir, DEFAULT_OUTPUT_FILE, DEFAULT_SERVER_URL);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Method description
     *
     * @param inputDir
     * @param outputFile
     * @param serverUrl
     * @throws Exception
     */
    public void process(File inputDir, String outputFile, String serverUrl) throws Exception {


        masterDocument = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();

        rootNode = masterDocument.createElement("Global");
        rootNode.setAttribute("name", "Data");
        rootNode.setAttribute("version", "1");
        masterDocument.appendChild(rootNode);
        if (epigenetics) {
            File[] files = inputDir.listFiles();
            Arrays.sort(files);
            doEpigenetics(files, serverUrl);
        } else {
            lookForChildren(inputDir, serverUrl);
        }

        // Transform document into XML
        TransformerFactory factory = TransformerFactory.newInstance();
        Transformer transformer = factory.newTransformer();
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

        StreamResult streamResult = new StreamResult(new StringWriter());
        DOMSource source = new DOMSource(masterDocument);
        transformer.transform(source, streamResult);


        // Create output dir
        String xmlString = streamResult.getWriter().toString();
        File tempXMLFile = new File(outputFile);
        FileWriter fileWriter = null;

        fileWriter = new FileWriter(tempXMLFile);
        fileWriter.write(xmlString);


        fileWriter.close();

    }

    private void createSampleInfoFiles(File[] files, String sampleInfoFile,
                                       String trackPropertiesFile) {
        Arrays.sort(files);
        try {
            PrintWriter sampleInfoWriter =
                    new PrintWriter(new BufferedWriter(new FileWriter(sampleInfoFile)));
            PrintWriter trackWriter =
                    new PrintWriter(new BufferedWriter(new FileWriter(trackPropertiesFile)));

            sampleInfoWriter.println("Name\tCell Type\tMark\t#height\t#renderer\t#color\t#min\t#max");

            for (File file : files) {
                String[] tokens = file.getName().split("\\.");
                if (tokens.length > 2) {
                    String trackName = file.getName().replace(".h5", "");
                    String cellType = tokens[0].replace("r1", "").replace("r2", "").replace("r3", "");
                    String mark = tokens[1];

                    sampleInfoWriter.print(trackName + "\t" + cellType + "\t" + mark);

                    int min = 0;
                    int max = 10;
                    String colorString = ColorUtilities.colorToString(Color.gray);
                    if (mark.toUpperCase().contains("K4")) {
                        max = 25;
                        colorString = "(0,150,0)";
                    } else if (mark.toUpperCase().contains("K9")) {
                        colorString = "(100,0,0)";
                    } else if (mark.toUpperCase().contains("K27")) {
                        colorString = "(255,0,0)";
                    }

                    sampleInfoWriter.println("\t40\tBAR_CHART\t" + colorString + "\t" + min
                            + "\t" + max);

                }
            }
            sampleInfoWriter.close();

        }
        catch (IOException iOException) {
        }
    }

    private void doEpigenetics(File[] files, String serverUrl) throws IOException {

        LinkedHashSet<String> cellTypes = new LinkedHashSet();
        LinkedHashSet<String> markers = new LinkedHashSet();
        Map<String, LinkedHashSet<File>> cellMarkerFileMap = new LinkedHashMap();

        for (File file : files) {
            String[] tokens = file.getName().split("\\.");

            if (tokens.length > 2) {
                String cellType = tokens[0];
                String marker = tokens[1];
                cellTypes.add(cellType);
                markers.add(marker);

                String cellMarker = cellType + "." + marker;
                LinkedHashSet<File> cellMarkerFiles = cellMarkerFileMap.get(cellMarker);
                if (cellMarkerFiles == null) {
                    cellMarkerFiles = new LinkedHashSet();
                    cellMarkerFileMap.put(cellMarker, cellMarkerFiles);
                }
                cellMarkerFiles.add(file);
            }
        }

        Element cellTypeRootNode = masterDocument.createElement("Category");
        cellTypeRootNode.setAttribute("name", "Cell Type");
        rootNode.appendChild(cellTypeRootNode);

        List<String> sortedCellTypes = new ArrayList(cellTypes);
        Collections.sort(sortedCellTypes);
        List<String> sortedMarks = new ArrayList(markers);
        Collections.sort(sortedMarks);

        for (String cellType : sortedCellTypes) {
            Element cellTypeNode = masterDocument.createElement("Category");
            cellTypeNode.setAttribute("name", cellType);
            cellTypeRootNode.appendChild(cellTypeNode);
            for (String marker : sortedMarks) {
                String cmKey = cellType + "." + marker;
                outputResources(cellMarkerFileMap, serverUrl, cellTypeNode, cmKey);
            }
        }

        Element markerRootNode = masterDocument.createElement("Category");
        markerRootNode.setAttribute("name", "Chromatin Mark");
        rootNode.appendChild(markerRootNode);
        for (String marker : sortedMarks) {
            Element markerNode = masterDocument.createElement("Category");
            markerNode.setAttribute("name", marker);
            markerRootNode.appendChild(markerNode);
            for (String cellType : sortedCellTypes) {
                String cmKey = cellType + "." + marker;
                outputResources(cellMarkerFileMap, serverUrl, markerNode, cmKey);
            }
        }

    }

    private void outputResources(Map<String, LinkedHashSet<File>> cellMarkerFileMap,
                                 String serverUrl, Element cellTypeRootNode, String cmKey)
            throws IOException, DOMException {
        LinkedHashSet<File> files = cellMarkerFileMap.get(cmKey);
        if (files != null) {
            System.out.println(cmKey + " -> " + files.size());
            if (files.size() > 0) {
                for (File file : files) {
                    String canonicalPath = file.getCanonicalPath();
                    Element resourceNode = masterDocument.createElement("Resource");
                    resourceNode.setAttribute("name", cmKey);
                    resourceNode.setAttribute("path", canonicalPath);
                    resourceNode.setAttribute("serverURL", serverUrl);
                    cellTypeRootNode.appendChild(resourceNode);
                }
            }
        }
    }

    private void lookForChildren(File dir, String serverUrl) throws IOException {

        if (dir.isDirectory()) {
            if (rootNode != null) {
                for (File file : dir.listFiles()) {
                    String canonicalPath = file.getCanonicalPath();
                    Element resourceNode = masterDocument.createElement("Resource");
                    resourceNode.setAttribute("name", file.getName());
                    resourceNode.setAttribute("path", canonicalPath);
                    resourceNode.setAttribute("serverURL", serverUrl);
                    rootNode.appendChild(resourceNode);
                }
            }
        } else {
            System.err.println("Input must be a directory: " + dir.getAbsolutePath());
        }
    }
}
