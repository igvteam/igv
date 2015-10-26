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

package org.broad.igv.session;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.session.IGVSessionReader.SessionAttribute;
import org.broad.igv.session.IGVSessionReader.SessionElement;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.TrackFilterElement;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import javax.swing.*;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.namespace.QName;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 */
public class SessionWriter {

    static Logger log = Logger.getLogger(SessionWriter.class);

    Session session;
    private static int CURRENT_VERSION = 8;

    private static final String TRACK_TAG = SessionElement.TRACK.getText();

    /**
     * Method description
     *
     * @param session
     * @param outputFile
     * @throws IOException
     */
    public void saveSession(Session session, File outputFile) throws IOException {

        if (session == null) {
            RuntimeException e = new RuntimeException("No session found to save!");
            log.error("Session Management Error", e);
        }

        this.session = session;

        if (outputFile == null) {
            RuntimeException e = new RuntimeException("Can't save session file: " + outputFile);
            log.error("Session Management Error", e);
        }

        String xmlString = createXmlFromSession(session, outputFile);

        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(outputFile);
            fileWriter.write(xmlString);
        } finally {
            if (fileWriter != null) {
                fileWriter.close();
            }
        }
    }


    public String createXmlFromSession(Session session, File outputFile) throws RuntimeException {

        String xmlString = null;

        try {

            // Create a DOM document
            DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document document = documentBuilder.newDocument();
            document.setStrictErrorChecking(true);

            // Global root element
            Element globalElement = document.createElement(SessionElement.SESSION.getText());

            globalElement.setAttribute(SessionAttribute.VERSION.getText(), String.valueOf(CURRENT_VERSION));

            String genomeId = GenomeManager.getInstance().getGenomeId();
            if (genomeId != null) {
                globalElement.setAttribute(SessionAttribute.GENOME.getText(), genomeId);
            }

            String locus = session.getLocusString();
            if (locus != null && !FrameManager.isGeneListMode()) {
                globalElement.setAttribute(SessionAttribute.LOCUS.getText(), locus);
            }

            String groupBy = IGV.getInstance().getGroupByAttribute();
            if (groupBy != null) {
                globalElement.setAttribute(SessionAttribute.GROUP_TRACKS_BY.getText(), groupBy);
            }

            int nextAutoscaleGroup = session.getNextAutoscaleGroup();
            if(nextAutoscaleGroup > 1) {
                globalElement.setAttribute(SessionAttribute.NEXT_AUTOSCALE_GROUP.getText(), String.valueOf(nextAutoscaleGroup));
            }

            if (session.isRemoveEmptyPanels()) {
                globalElement.setAttribute("removeEmptyTracks", "true");
            }

            globalElement.setAttribute(SessionAttribute.HAS_GENE_TRACK.getText(), "" + IGV.getInstance().hasGeneTrack());
            globalElement.setAttribute(SessionAttribute.HAS_SEQ_TRACK.getText(), "" + IGV.getInstance().hasSequenceTrack());

            globalElement.setAttribute("path", outputFile.getAbsolutePath());


            // Resource Files
            writeResources(outputFile, globalElement, document);

            // Panels
            writePanels(globalElement, document);

            // Panel layout
            writePanelLayout(globalElement, document);


            // Regions of Interest
            writeRegionsOfInterest(globalElement, document);

            // Filter
            writeFilters(session, globalElement, document);

            if (FrameManager.isGeneListMode()) {
                writeGeneList(globalElement, document);
            }

            // Hidden attributes
            if (session.getHiddenAttributes() != null && session.getHiddenAttributes().size() > 0) {
                writeHiddenAttributes(session, globalElement, document);
            }


            document.appendChild(globalElement);

            // Transform document into XML
            xmlString = Utilities.getString(document);
        } catch (Exception e) {
            String message = "An error has occurred while trying to create the session!";
            log.error(message, e);
            JOptionPane.showMessageDialog(IGV.getMainFrame(), message);
            throw new RuntimeException(e);
        }

        return xmlString;
    }


    private void writeFilters(Session session, Element globalElement, Document document) {
        TrackFilter trackFilter = session.getFilter();
        if (trackFilter != null) {

            Element filter = document.createElement(SessionElement.FILTER.getText());

            filter.setAttribute(SessionAttribute.NAME.getText(), trackFilter.getName());

            if (IGV.getInstance().isFilterMatchAll()) {
                filter.setAttribute(SessionAttribute.FILTER_MATCH.getText(), "all");
            } else if (!IGV.getInstance().isFilterMatchAll()) {
                filter.setAttribute(SessionAttribute.FILTER_MATCH.getText(), "any");
            } else {    // Defaults to match all
                filter.setAttribute(SessionAttribute.FILTER_MATCH.getText(), "all");
            }

            if (IGV.getInstance().isFilterShowAllTracks()) {
                filter.setAttribute(SessionAttribute.FILTER_SHOW_ALL_TRACKS.getText(), "true");
            } else {    // Defaults
                filter.setAttribute(SessionAttribute.FILTER_SHOW_ALL_TRACKS.getText(), "false");
            }
            globalElement.appendChild(filter);

            // Process FilterElement elements
            Iterator iterator = session.getFilter().getFilterElements();
            while (iterator.hasNext()) {

                TrackFilterElement trackFilterElement = (TrackFilterElement) iterator.next();

                Element filterElementElement =
                        document.createElement(SessionElement.FILTER_ELEMENT.getText());
                filterElementElement.setAttribute(SessionAttribute.ITEM.getText(),
                        trackFilterElement.getSelectedItem());
                filterElementElement.setAttribute(
                        SessionAttribute.OPERATOR.getText(),
                        trackFilterElement.getComparisonOperator().getValue());
                filterElementElement.setAttribute(SessionAttribute.VALUE.getText(),
                        trackFilterElement.getValue());
                filterElementElement.setAttribute(
                        SessionAttribute.BOOLEAN_OPERATOR.getText(),
                        trackFilterElement.getBooleanOperator().getValue());
                filter.appendChild(filterElementElement);
            }
        }
    }

    private void writeRegionsOfInterest(Element globalElement, Document document) {
        Collection<RegionOfInterest> regions = session.getAllRegionsOfInterest();
        if ((regions != null) && !regions.isEmpty()) {

            Element regionsElement = document.createElement(SessionElement.REGIONS.getText());
            for (RegionOfInterest region : regions) {
                Element regionElement = document.createElement(SessionElement.REGION.getText());
                regionElement.setAttribute(SessionAttribute.CHROMOSOME.getText(), region.getChr());
                regionElement.setAttribute(SessionAttribute.START_INDEX.getText(), String.valueOf(region.getStart()));
                regionElement.setAttribute(SessionAttribute.END_INDEX.getText(), String.valueOf(region.getEnd()));
                if (region.getDescription() != null) {
                    regionElement.setAttribute(SessionAttribute.DESCRIPTION.getText(), region.getDescription());
                }
                regionsElement.appendChild(regionElement);
            }
            globalElement.appendChild(regionsElement);
        }
    }

    private void writeHiddenAttributes(Session session, Element globalElement, Document document) {
        Element hiddenAttributes = document.createElement(SessionElement.HIDDEN_ATTRIBUTES.getText());
        for (String attribute : session.getHiddenAttributes()) {
            Element regionElement = document.createElement(SessionElement.ATTRIBUTE.getText());
            regionElement.setAttribute(IGVSessionReader.SessionAttribute.NAME.getText(), attribute);
            hiddenAttributes.appendChild(regionElement);
        }
        globalElement.appendChild(hiddenAttributes);

    }

    private void writeGeneList(Element globalElement, Document document) {

        GeneList geneList = session.getCurrentGeneList();

        if (geneList != null) {

            Element geneListElement = document.createElement(SessionElement.GENE_LIST.getText());
            geneListElement.setAttribute(IGVSessionReader.SessionAttribute.NAME.getText(), geneList.getName());

            StringBuffer genes = new StringBuffer();
            for (String gene : geneList.getLoci()) {
                genes.append(gene);
                genes.append("\n");
            }

            geneListElement.setTextContent(genes.toString());

            globalElement.appendChild(geneListElement);


            // Now store the list of frames visible
            for (ReferenceFrame frame : FrameManager.getFrames()) {

                Element frameElement = document.createElement(SessionElement.FRAME.getText());
                frameElement.setAttribute(IGVSessionReader.SessionAttribute.NAME.getText(), frame.getName());
                frameElement.setAttribute(IGVSessionReader.SessionAttribute.CHR.getText(), frame.getChrName());
                frameElement.setAttribute(IGVSessionReader.SessionAttribute.START.getText(), String.valueOf(frame.getOrigin()));
                frameElement.setAttribute(IGVSessionReader.SessionAttribute.END.getText(), String.valueOf(frame.getEnd()));

                geneListElement.appendChild(frameElement);

            }
        }
    }

    private void writeResources(File outputFile, Element globalElement, Document document) throws IOException {

        Collection<ResourceLocator> resourceLocators = getResourceLocatorSet();

        if ((resourceLocators != null) && !resourceLocators.isEmpty()) {

            Element filesElement = document.createElement(SessionElement.RESOURCES.getText());

            for (ResourceLocator resourceLocator : resourceLocators) {
                if (resourceLocator.exists() || !(resourceLocator.getPath() == null)) {

                    //RESOURCE ELEMENT
                    Element dataFileElement = document.createElement(SessionElement.RESOURCE.getText());

                    //REQUIRED ATTRIBUTES - Cannot be null

                    boolean useRelative = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SESSION_RELATIVE_PATH);

                    String relativePath = useRelative ?
                            FileUtils.getRelativePath(outputFile.getAbsolutePath(), resourceLocator.getPath()) :
                            resourceLocator.getPath();

                    dataFileElement.setAttribute(SessionAttribute.PATH.getText(), relativePath);

                    //OPTIONAL ATTRIBUTES

                    if (resourceLocator.getName() != null) {
                        dataFileElement.setAttribute(SessionAttribute.NAME.getText(), resourceLocator.getName());
                    }
                    if (resourceLocator.getDBUrl() != null) {
                        dataFileElement.setAttribute(SessionAttribute.SERVER_URL.getText(), resourceLocator.getDBUrl());
                    }
                    if (resourceLocator.getTrackInfoURL() != null) {
                        dataFileElement.setAttribute(SessionAttribute.HYPERLINK.getText(), resourceLocator.getTrackInfoURL());
                    }
                    if (resourceLocator.getFeatureInfoURL() != null) {
                        dataFileElement.setAttribute(SessionAttribute.FEATURE_URL.getText(), resourceLocator.getFeatureInfoURL());
                    }
                    if (resourceLocator.getDescription() != null) {
                        dataFileElement.setAttribute(SessionAttribute.DESCRIPTION.getText(), resourceLocator.getDescription());
                    }
                    if (resourceLocator.getType() != null) {
                        dataFileElement.setAttribute(SessionAttribute.TYPE.getText(), resourceLocator.getType());
                    }
                    if (resourceLocator.getIndexPath() != null) {
                        dataFileElement.setAttribute(SessionAttribute.INDEX.getText(), resourceLocator.getIndexPath());
                    }
                    if (resourceLocator.getCoverage() != null) {
                        dataFileElement.setAttribute(SessionAttribute.COVERAGE.getText(), resourceLocator.getCoverage());
                    }
                    if (resourceLocator.getTrackLine() != null) {
                        dataFileElement.setAttribute(SessionAttribute.TRACK_LINE.getText(), resourceLocator.getTrackLine());
                    }
                    filesElement.appendChild(dataFileElement);
                }
            }
            globalElement.appendChild(filesElement);
        }
    }

    private void writePanels(Element globalElement, Document document) throws DOMException {

        for (TrackPanel trackPanel : IGV.getInstance().getTrackPanels()) {

            // TODO -- loop through panels groups, rather than skipping groups to tracks

            List<Track> tracks = trackPanel.getTracks();
            if ((tracks != null) && !tracks.isEmpty()) {

                Element panelElement = document.createElement(SessionElement.PANEL.getText());
                panelElement.setAttribute("name", trackPanel.getName());
                panelElement.setAttribute("height", String.valueOf(trackPanel.getHeight()));
                panelElement.setAttribute("width", String.valueOf(trackPanel.getWidth()));

                //We create a temporary element into which to marshall, so we
                //can add custom attributes
                Element tmpTrackParent = document.createElement("dummy");

                try {
                    Marshaller m = IGVSessionReader.getJAXBContext().createMarshaller();
                    m.setProperty(Marshaller.JAXB_FRAGMENT, true);
                    for (Track track : tracks) {

                        marshalTrack(m, track, tmpTrackParent, track.getClass());

                        Element trackElement = (Element) tmpTrackParent.getChildNodes().item(0);

                        for (Map.Entry<String, String> attrValue : track.getPersistentState().entrySet()) {
                            trackElement.setAttribute(attrValue.getKey(), attrValue.getValue());
                        }

                        marshalTrackChildren(m, track, trackElement);

                        panelElement.appendChild(trackElement);
                    }
                } catch (JAXBException e) {
                    throw new RuntimeException(e);
                }
                globalElement.appendChild(panelElement);
            }
        }
    }

    /**
     * Attempt to marshall the {@code track} into {@code trackParent} as it's
     * own class, if that fails, try the superclass, and so on up
     *
     * @param m
     * @param track
     * @param trackParent
     * @throws javax.xml.bind.JAXBException
     */
    private static void marshalTrack(Marshaller m, Track track, Node trackParent, Class marshalClass) throws JAXBException {

        if (marshalClass == null || marshalClass.equals(Object.class)) {
            throw new JAXBException(track.getClass() + " and none of its superclasses are known");
        }

        if (AbstractTrack.knownUnknownTrackClasses.contains(marshalClass)) {
            marshalTrack(m, track, trackParent, marshalClass.getSuperclass());
            return;
        }

        JAXBElement el;
        try {
            el = new JAXBElement(new QName("", TRACK_TAG), marshalClass, track);
            m.marshal(el, trackParent);
        } catch (JAXBException e) {
            AbstractTrack.knownUnknownTrackClasses.add(marshalClass);
            marshalTrack(m, track, trackParent, marshalClass.getSuperclass());
        }
    }


    /**
     * Because we are using JAXB piecewise, and also because JAXB can't handle
     * interfaces, we marshal certain track children here
     *
     * @param m
     * @param track
     * @param trackElement
     */
    private void marshalTrackChildren(Marshaller m, Track track, Element trackElement) throws JAXBException {
        if (track instanceof FeatureTrack) {
            FeatureTrack featureTrack = (FeatureTrack) track;
            featureTrack.marshalSource(m, trackElement);
        } else if (track instanceof DataSourceTrack) {
            DataSourceTrack dataSourceTrack = (DataSourceTrack) track;
            dataSourceTrack.marshalSource(m, trackElement);
        }
    }


    private void writePanelLayout(Element globalElement, Document document) {

        double[] dividerFractions = IGV.getInstance().getMainPanel().getDividerFractions();
        if (dividerFractions.length > 0) {

            Element panelLayout = document.createElement(SessionElement.PANEL_LAYOUT.getText());
            globalElement.appendChild(panelLayout);

            StringBuffer locString = new StringBuffer();
            locString.append(String.valueOf(dividerFractions[0]));
            for (int i = 1; i < dividerFractions.length; i++) {
                locString.append("," + dividerFractions[i]);
            }
            panelLayout.setAttribute("dividerFractions", locString.toString());


        }

    }

    /**
     * @return A set of the load data files.
     */
    public Collection<ResourceLocator> getResourceLocatorSet() {

        Collection<ResourceLocator> locators = new ArrayList();

        Collection<ResourceLocator> currentTrackFileLocators =
                IGV.getInstance().getDataResourceLocators();

        if (currentTrackFileLocators != null) {
            for (ResourceLocator locator : currentTrackFileLocators) {
                locators.add(locator);
            }
        }

        Collection<ResourceLocator> loadedAttributeResources =
                AttributeManager.getInstance().getLoadedResources();

        if (loadedAttributeResources != null) {
            for (ResourceLocator attributeLocator : loadedAttributeResources) {
                locators.add(attributeLocator);
            }
        }

        return locators;
    }

}

