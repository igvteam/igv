/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.hic.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.hic.MainWindow;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.swing.*;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 5/8/12
 */
public class HiCTrackManager {

    static String path = "http://www.broadinstitute.org/igvdata/hic/tracksMenu.xml";
    //static String path = "/Users/jrobinso/Documents/IGV/hg19_encode.xml";

    // Category => list of locators
    private static Map<String, List<ResourceLocator>> categoryLocatorMap = null;

    // track name => locator
    private static Map<String, ResourceLocator> locatorMap;

    private static java.util.List<Track> loadedTracks = new ArrayList();

    public static void loadTrackFromFile(final MainWindow parent) {
        FileDialog dlg = new FileDialog(parent);
        dlg.setMode(FileDialog.LOAD);
        dlg.setVisible(true);
        String file = dlg.getFile();
        if (file != null) {

            final File f = new File(dlg.getDirectory(), dlg.getFile());
            // Load new tracks
            Runnable runnable = new Runnable() {
                public void run() {

                    Genome genome = loadGenome();

                    String path = f.getAbsolutePath();
                    // genome = GenomeManager.getInstance().getCurrentGenome();
                    ResourceLocator locator = new ResourceLocator(path);
                    List<Track> tracks = (new TrackLoader()).load(locator, genome);
                    loadedTracks.addAll(tracks);
                    parent.updateTrackPanel();

                }
            };
            parent.executeLongRunningTask(runnable);
        }
    }

    public static void loadHostedTrack(final MainWindow parent) {

        Genome genome = loadGenome();

        Map<String, List<ResourceLocator>> locators = getTrackLocators();
        HiCLoadDialog dlg = new HiCLoadDialog(parent, locators, loadedTracks);
        dlg.setVisible(true);

        if (!dlg.isCanceled()) {
            final Collection<String> selectedTracks = dlg.getSelectedTracks();

            // Unload de-selected tracks
            boolean trackRemoved = false;
            Iterator<Track> trackIterator = loadedTracks.iterator();
            final Set<String> loadedTrackNames = new HashSet<String>();
            while (trackIterator.hasNext()) {
                final Track track = trackIterator.next();
                if (!selectedTracks.contains(track.getName())) {
                    trackIterator.remove();
                    trackRemoved = true;
                } else {
                    loadedTrackNames.add(track.getName());
                }
            }

            if (trackRemoved) {
                parent.updateTrackPanel();
            }


            // Load new tracks
            if (selectedTracks.size() > 0) {
                Runnable runnable = new Runnable() {
                    public void run() {
                        for (String trackName : selectedTracks) {
                            if (!loadedTrackNames.contains(trackName)) {
                                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                                ResourceLocator locator = locatorMap.get(trackName);
                                List<Track> tracks = (new TrackLoader()).load(locator, genome);

                                // If the track supports a 90% window function use it.
                                for(Track t : tracks) {
                                    Collection<WindowFunction> wfs =  t.getAvailableWindowFunctions();
                                    if(wfs != null && wfs.contains(WindowFunction.percentile90)) {
                                        t.setWindowFunction(WindowFunction.percentile90);
                                    }
                                }

                                loadedTracks.addAll(tracks);
                            }
                        }
                        parent.updateTrackPanel();

                    }
                };
                parent.executeLongRunningTask(runnable);
            }
        }

    }


    private static Genome loadGenome() {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome == null) {
            String genomePath = "http://igvdata.broadinstitute.org/genomes/hg19.genome";
            try {
                genome = GenomeManager.getInstance().loadGenome(genomePath, null);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return genome;
    }

    public static List<Track> getLoadedTracks() {
        return loadedTracks;
    }


    public static synchronized Map<String, List<ResourceLocator>> getTrackLocators() {
        if (categoryLocatorMap == null) {
            try {
                categoryLocatorMap = parseResourceFile(path);
                locatorMap = new HashMap<String, ResourceLocator>();
                for (java.util.List<ResourceLocator> locatorList : categoryLocatorMap.values()) {
                    for (ResourceLocator loc : locatorList) {
                        locatorMap.put(loc.getTrackName(), loc);
                    }
                }

            } catch (ParserConfigurationException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (SAXException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        return categoryLocatorMap;
    }

    private static Map<String, List<ResourceLocator>> parseResourceFile(String file) throws ParserConfigurationException, IOException, SAXException {

        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        Document doc = dBuilder.parse(file);

        Map<String, List<ResourceLocator>> locators = Collections.synchronizedMap(new LinkedHashMap<String, List<ResourceLocator>>());

        NodeList categoryNodes = doc.getElementsByTagName("Category");
        int nNodes = categoryNodes.getLength();
        for (int i = 0; i < nNodes; i++) {
            Node node = categoryNodes.item(i);
            String nodeName = node.getAttributes().getNamedItem("name").getNodeValue();

            NodeList resourceNodes = node.getChildNodes();
            List<ResourceLocator> locatorList = new ArrayList(resourceNodes.getLength());
            for (int j = 0; j < resourceNodes.getLength(); j++) {
                Node resourceNode = resourceNodes.item(j);
                if (resourceNode.getNodeName().equals("Resource")) {
                    NamedNodeMap attributes = resourceNode.getAttributes();
                    Node nameNode = attributes.getNamedItem("name");
                    Node pathNode = attributes.getNamedItem("path");
                    Node trackLineNode = attributes.getNamedItem("trackLine");
                    Node colorNode = attributes.getNamedItem("color");
                    if (nameNode == null || pathNode == null) {
                        System.out.println("Skipping " + node.toString());
                        continue;
                    }

                    ResourceLocator rl = new ResourceLocator(pathNode.getNodeValue());
                    rl.setName(nameNode.getNodeValue());

                    if (trackLineNode != null) {
                        rl.setTrackLine(trackLineNode.getNodeValue());
                    }
                    if (colorNode != null) {
                        try {
                            rl.setColor(ColorUtilities.stringToColor(colorNode.getNodeValue()));
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }

                    locatorList.add(rl);

                }
            }

            locators.put(nodeName, locatorList);
        }
        return locators;


    }


}
