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
package org.broad.igv.session;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.TrackFilterElement;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.FilterElement.BooleanOperator;
import org.broad.igv.util.FilterElement.Operator;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.List;

/**
 *
 */
public class SessionReader {

    private static Logger log = Logger.getLogger(SessionReader.class);
    private static String INPUT_FILE_KEY = "INPUT_FILE_KEY";
    // Temporary values used in processing

    private Collection<ResourceLocator> dataFiles;
    private Collection<ResourceLocator> missingDataFiles;
    private static Map<String, String> attributeSynonymMap = new HashMap();
    private boolean panelElementPresent = false;
    private int version;

    private IGV igv;


    /**
     * Map of track id -> track.  It is important to maintin the order in which tracks are added, thus
     * the use of LinkedHashMap.
     */
    Map<String, List<Track>> trackDictionary = new LinkedHashMap();
    private Track geneTrack = null;
    private Track seqTrack = null;


    static {
        attributeSynonymMap.put("DATA FILE", "DATA SET");
        attributeSynonymMap.put("TRACK NAME", "NAME");
    }

    /**
     * Session Element types
     */
    public static enum SessionElement {


        PANEL("Panel"),
        PANEL_LAYOUT("PanelLayout"),
        TRACK("Track"),
        COLOR_SCALE("ColorScale"),
        COLOR_SCALES("ColorScales"),
        DATA_TRACK("DataTrack"),
        DATA_TRACKS("DataTracks"),
        FEATURE_TRACKS("FeatureTracks"),
        DATA_FILE("DataFile"),
        RESOURCE("Resource"),
        RESOURCES("Resources"),
        FILES("Files"),
        FILTER_ELEMENT("FilterElement"),
        FILTER("Filter"),
        SESSION("Session"),
        GLOBAL("Global"),
        REGION("Region"),
        REGIONS("Regions"),
        DATA_RANGE("DataRange"),
        PREFERENCES("Preferences"),
        PROPERTY("Property"),
        GENE_LIST("GeneList"),
        HIDDEN_ATTRIBUTES("HiddenAttributes"),
        ATTRIBUTE("Attribute"),
        FRAME("Frame");

        private String name;

        SessionElement(String name) {
            this.name = name;
        }

        public String getText() {
            return name;
        }

        @Override
        public String toString() {
            return getText();
        }

        static public SessionElement findEnum(String value) {

            if (value == null) {
                return null;
            } else {
                return SessionElement.valueOf(value);
            }
        }
    }

    /**
     * Session Attribute types
     */
    public static enum SessionAttribute {

        BOOLEAN_OPERATOR("booleanOperator"),
        COLOR("color"),
        ALT_COLOR("altColor"),
        COLOR_MODE("colorMode"),
        CHROMOSOME("chromosome"),
        END_INDEX("end"),
        EXPAND("expand"),
        SQUISH("squish"),
        DISPLAY_MODE("displayMode"),
        FILTER_MATCH("match"),
        FILTER_SHOW_ALL_TRACKS("showTracks"),
        GENOME("genome"),
        GROUP_TRACKS_BY("groupTracksBy"),
        HEIGHT("height"),
        ID("id"),
        ITEM("item"),
        LOCUS("locus"),
        NAME("name"),
        SAMPLE_ID("sampleID"),
        RESOURCE_TYPE("resourceType"),
        OPERATOR("operator"),
        RELATIVE_PATH("relativePath"),
        RENDERER("renderer"),
        SCALE("scale"),
        START_INDEX("start"),
        VALUE("value"),
        VERSION("version"),
        VISIBLE("visible"),
        WINDOW_FUNCTION("windowFunction"),
        RENDER_NAME("renderName"),
        GENOTYPE_HEIGHT("genotypeHeight"),
        VARIANT_HEIGHT("variantHeight"),
        PREVIOUS_HEIGHT("previousHeight"),
        FEATURE_WINDOW("featureVisibilityWindow"),
        DISPLAY_NAME("displayName"),
        COLOR_SCALE("colorScale"),

        //RESOURCE ATTRIBUTES
        PATH("path"),
        LABEL("label"),
        SERVER_URL("serverURL"),
        HYPERLINK("hyperlink"),
        INFOLINK("infolink"),
        URL("url"),
        FEATURE_URL("featureURL"),
        DESCRIPTION("description"),
        TYPE("type"),
        COVERAGE("coverage"),
        TRACK_LINE("trackLine"),

        CHR("chr"),
        START("start"),
        END("end");

        //TODO Add the following into the Attributes
        /*
        boolean shadeBases;
        boolean shadeCenters;
        boolean flagUnmappedPairs;
        boolean showAllBases;
        int insertSizeThreshold;
        boolean colorByStrand;
        boolean colorByAmpliconStrand;
         */


        private String name;

        SessionAttribute(String name) {
            this.name = name;
        }

        public String getText() {
            return name;
        }

        @Override
        public String toString() {
            return getText();
        }

        static public SessionAttribute findEnum(String value) {

            if (value == null) {
                return null;
            } else {
                return SessionAttribute.valueOf(value);
            }
        }
    }


    public SessionReader(IGV igv) {
        this.igv = igv;
    }


    /**
     * @param inputStream
     * @param session
     *@param sessionName  @return
     * @throws RuntimeException
     */

    public void loadSession(InputStream inputStream, Session session, String sessionName)
            throws RuntimeException {


        log.debug("Load session");


        Document document = null;
        try {
            document = createDOMDocumentFromXmlFile(inputStream);
        } catch (Exception e) {
            log.error("Session Management Error", e);
            throw new RuntimeException(e);
        }

        HashMap additionalInformation = new HashMap();
        additionalInformation.put(INPUT_FILE_KEY, sessionName);

        NodeList nodes = document.getElementsByTagName(SessionElement.GLOBAL.getText());
        if (nodes == null || nodes.getLength() == 0) {
            nodes = document.getElementsByTagName(SessionElement.SESSION.getText());
        }

        processRootNode(session, nodes.item(0), additionalInformation);

        // Add tracks not explicitly set in file.  It is legal to define sessions with the DataFile section only (no
        // Panel or Track elements).
        addLeftoverTracks(trackDictionary.values());

        if (session.getGroupTracksBy() != null && session.getGroupTracksBy().length() > 0) {
            IGV.getInstance().getTrackManager().setGroupByAttribute(session.getGroupTracksBy());
        }

        if (session.isRemoveEmptyTracks()) {
            IGV.getInstance().getMainPanel().removeEmptyDataPanels();
        }

        IGV.getInstance().getTrackManager().resetOverlayTracks();

    }


    private void processRootNode(Session session, Node node, HashMap additionalInformation) {

        if ((node == null) || (session == null)) {
            MessageUtils.showMessage("Invalid session file: root node not found");
            return;
        }

        String nodeName = node.getNodeName();
        if (!(nodeName.equalsIgnoreCase(SessionElement.GLOBAL.getText()) || nodeName.equalsIgnoreCase(SessionElement.SESSION.getText()))) {
            MessageUtils.showMessage("Session files must begin with a \"Global\" or \"Session\" element.  Found: " + nodeName);
        }
        process(session, node, additionalInformation);

        Element element = (Element) node;
        session.setGenome(getAttribute(element, SessionAttribute.GENOME.getText()));
        session.setLocus(getAttribute(element, SessionAttribute.LOCUS.getText()));
        session.setGroupTracksBy(getAttribute(element, SessionAttribute.GROUP_TRACKS_BY.getText()));

        String removeEmptyTracks = getAttribute(element, "removeEmptyTracks");
        if (removeEmptyTracks != null) {
            try {
                Boolean b = Boolean.parseBoolean(removeEmptyTracks);
                session.setRemoveEmptyTracks(b);
            }
            catch (Exception e) {
                log.error("Error parsing removeEmptyTracks string: " + removeEmptyTracks, e);
            }
        }

        String versionString = getAttribute(element, SessionAttribute.VERSION.getText());
        try {
            version = Integer.parseInt(versionString);
        }
        catch (NumberFormatException e) {
            log.error("Non integer version number in session file: " + versionString);
        }
        session.setVersion(version);

        geneTrack = IGV.getInstance().getTrackManager().getGeneTrack();
        if (geneTrack != null) {
            trackDictionary.put(geneTrack.getId(), Arrays.asList(geneTrack));
        }
        seqTrack = IGV.getInstance().getTrackManager().getSequenceTrack();
        if (seqTrack != null) {
            trackDictionary.put(seqTrack.getId(), Arrays.asList(seqTrack));
        }


        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);

        // ReferenceFrame.getInstance().invalidateLocationScale();
    }

    //TODO Check to make sure tracks are not being created twice
    //TODO -- DONT DO THIS FOR NEW SESSIONS

    private void addLeftoverTracks(Collection<List<Track>> tmp) {
        if (version < 3 || !panelElementPresent) {
            for (List<Track> tracks : tmp) {
                for (Track track : tracks) {
                    if (track != geneTrack && track != seqTrack && track.getResourceLocator() != null) {
                        TrackPanel group = IGV.getInstance().getTrackManager().getPanelFor(track.getResourceLocator());
                        group.addTrack(track);
                    }
                }
            }
        }

    }


    /**
     * Process a single session element node.
     *
     * @param session
     * @param element
     */
    private void process(Session session, Node element, HashMap additionalInformation) {

        if ((element == null) || (session == null)) {
            return;
        }

        String nodeName = element.getNodeName();

        if (nodeName.equalsIgnoreCase(SessionElement.FILES.getText())) {
            processFiles(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.DATA_FILE.getText())) {
            processDataFile(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.RESOURCES.getText())) {
            processResources(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.RESOURCE.getText())) {
            processResource(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.REGIONS.getText())) {
            processRegions(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.REGION.getText())) {
            processRegion(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.GENE_LIST.getText())) {
            processGeneList(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.FILTER.getText())) {
            processFilter(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.FILTER_ELEMENT.getText())) {
            processFilterElement(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.COLOR_SCALES.getText())) {
            processColorScales(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.COLOR_SCALE.getText())) {
            processColorScale(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.PREFERENCES.getText())) {
            processPreferences(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.DATA_TRACKS.getText()) ||
                nodeName.equalsIgnoreCase(SessionElement.FEATURE_TRACKS.getText()) ||
                nodeName.equalsIgnoreCase(SessionElement.PANEL.getText())) {
            processPanel(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.PANEL_LAYOUT.getText())) {
            processPanelLayout(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.HIDDEN_ATTRIBUTES.getText())) {
            processHiddenAttributes(session, (Element) element, additionalInformation);
        }

    }


    /**
     * Process the Files element.
     * <p/>
     * The RELATIVE_PATH attribute specifies whether file paths are relative
     * or absolute.
     *
     * @param session
     * @param element
     */
    private void processFiles(Session session, Element element, HashMap additionalInformation) {

        dataFiles = new ArrayList();
        missingDataFiles = new ArrayList();
        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);

        if (missingDataFiles.size() > 0) {
            StringBuffer message = new StringBuffer();
            message.append("<html>The following data file(s) could not be located.<ul>");
            for (ResourceLocator file : missingDataFiles) {
                if (file.isLocal()) {
                    message.append("<li>");
                    message.append(file.getPath());
                    message.append("</li>");
                } else {
                    message.append("<li>Server: ");
                    message.append(file.getServerURL());
                    message.append("  Path: ");
                    message.append(file.getPath());
                    message.append("</li>");
                }
            }
            message.append("</ul>");
            message.append("Common reasons for this include: ");
            message.append("<ul><li>The session or data files have been moved.</li> ");
            message.append("<li>The data files are located on a drive that is not currently accessible.</li></ul>");
            message.append("</html>");

            MessageUtils.showMessage(message.toString());
        }
        if (dataFiles.size() > 0) {

            for (ResourceLocator locator : dataFiles) {

                List<Track> tracks = igv.getTrackManager().load(locator);

                for (Track track : tracks) {

                    String id = track.getId();
                    List<Track> trackList = trackDictionary.get(id);
                    if (trackList == null) {
                        trackList = new ArrayList();
                        trackDictionary.put(id, trackList);

                    }
                    trackList.add(track);
                }

            }
        }
        dataFiles = null;
    }

    /**
     * //TODO -- I think this method can be removed
     *
     * @param session
     * @param element
     * @Deprecated
     * @Deprectated -- user processResource
     * Process the data file element.  If relativePaths == true treat the
     * file path as relative to the session file path.  If false
     */

    private void processDataFile(Session session, Element element, HashMap additionalInformation) {

        ResourceLocator resourceLocator = null;
        String serverURL = getAttribute(element, SessionAttribute.SERVER_URL.getText());
        String filePath = getAttribute(element, SessionAttribute.NAME.getText());

        // e.g. DAS
        String resourceType = getAttribute(element, SessionAttribute.RESOURCE_TYPE.getText());

        // If file is local
        if ((serverURL == null || serverURL.trim().equals("")) &&
                !(IGVHttpUtils.isURL(filePath.toLowerCase()))) {
            String relPathValue = getAttribute(element, SessionAttribute.RELATIVE_PATH.getText());
            boolean relativePaths = ((relPathValue != null) && relPathValue.equalsIgnoreCase("true"));
            File parent = (relativePaths ? new File(session.getPath()).getParentFile() : null);
            File file = new File(parent, filePath);
            resourceLocator = new ResourceLocator(file.getAbsolutePath());
        } else {    // ...else must be from Server
            resourceLocator = new ResourceLocator(serverURL, filePath);
        }

        if (resourceType != null) {
            resourceLocator.setType(resourceType);
        }

        if (resourceLocator.exists()) {
            dataFiles.add(resourceLocator);
        } else {
            missingDataFiles.add(resourceLocator);
        }

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);
    }

    private void processResources(Session session, Element element, HashMap additionalInformation) {
        processFiles(session, element, additionalInformation);
    }

    private void processResource(Session session, Element element, HashMap additionalInformation) {

        String label = getAttribute(element, SessionAttribute.LABEL.getText());
        String name = getAttribute(element, SessionAttribute.NAME.getText());
        String sampleId = getAttribute(element, SessionAttribute.SAMPLE_ID.getText());
        String description = getAttribute(element, SessionAttribute.DESCRIPTION.getText());
        String type = getAttribute(element, SessionAttribute.TYPE.getText());
        String coverage = getAttribute(element, SessionAttribute.COVERAGE.getText());
        String trackLine = getAttribute(element, SessionAttribute.TRACK_LINE.getText());
        String colorString = getAttribute(element, SessionAttribute.COLOR.getText());

        String relPathValue = getAttribute(element, SessionAttribute.RELATIVE_PATH.getText());
        boolean relativePaths = ((relPathValue != null) && relPathValue.equalsIgnoreCase("true"));
        String serverURL = getAttribute(element, SessionAttribute.SERVER_URL.getText());
        String path = getAttribute(element, SessionAttribute.PATH.getText());

        ResourceLocator resourceLocator = new ResourceLocator(serverURL, path);
        if (relativePaths) {
            if (FileUtils.isRemote(session.getPath())) {
                int idx = session.getPath().lastIndexOf("/");
                String basePath = session.getPath().substring(0, idx);
                String resPath = basePath + "/" + path;
                resourceLocator = new ResourceLocator(serverURL, resPath);
            } else {
                File parent = (relativePaths ? new File(session.getPath()).getParentFile() : null);
                File file = new File(parent, path);
                resourceLocator = new ResourceLocator(serverURL, file.getAbsolutePath());
            }
        }


        String url = getAttribute(element, SessionAttribute.URL.getText());
        if (url == null) {
            url = getAttribute(element, SessionAttribute.FEATURE_URL.getText());
        }
        resourceLocator.setUrl(url);

        String infolink = getAttribute(element, SessionAttribute.HYPERLINK.getText());
        if (infolink == null) {
            infolink = getAttribute(element, SessionAttribute.INFOLINK.getText());
        }
        resourceLocator.setInfolink(infolink);


        // Label is deprecated in favor of name.
        if (name != null) {
            resourceLocator.setName(name);
        } else {
            resourceLocator.setName(label);
        }

        resourceLocator.setSampleId(sampleId);


        resourceLocator.setDescription(description);
        // This test added to get around earlier bug in the writer
        if (type != null && !type.equals("local")) {
            resourceLocator.setType(type);
        }
        resourceLocator.setCoverage(coverage);
        resourceLocator.setTrackLine(trackLine);

        if (colorString != null) {
            try {
                Color c = ColorUtilities.stringToColor(colorString);
                resourceLocator.setColor(c);
            } catch (Exception e) {
                log.error("Error setting color: ", e);
            }
        }

        if (resourceLocator.exists()) {
            dataFiles.add(resourceLocator);
        } else {
            missingDataFiles.add(resourceLocator);
        }

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);

    }

    private void processRegions(Session session, Element element, HashMap additionalInformation) {

        session.clearRegionsOfInterest();
        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);
    }

    private void processRegion(Session session, Element element, HashMap additionalInformation) {

        String chromosome = getAttribute(element, SessionAttribute.CHROMOSOME.getText());
        String start = getAttribute(element, SessionAttribute.START_INDEX.getText());
        String end = getAttribute(element, SessionAttribute.END_INDEX.getText());
        String description = getAttribute(element, SessionAttribute.DESCRIPTION.getText());

        RegionOfInterest region = new RegionOfInterest(chromosome, new Integer(start), new Integer(end), description);
        IGV.getInstance().addRegionOfInterest(region);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);
    }


    private void processHiddenAttributes(Session session, Element element, HashMap additionalInformation) {

//        session.clearRegionsOfInterest();
        NodeList elements = element.getChildNodes();
        if (elements.getLength() > 0) {
            Set<String> attributes = new HashSet();
            for (int i = 0; i < elements.getLength(); i++) {
                Node childNode = elements.item(i);
                if (childNode.getNodeName().equals(SessionReader.SessionElement.ATTRIBUTE.getText())) {
                    attributes.add(((Element) childNode).getAttribute(SessionReader.SessionAttribute.NAME.getText()));
                }
            }
            session.setHiddenAttributes(attributes);
        }
    }

    private void processGeneList(Session session, Element element, HashMap additionalInformation) {

        String name = getAttribute(element, SessionAttribute.NAME.getText());

        String txt = element.getTextContent();
        String[] genes = txt.trim().split("\\s+");
        GeneList gl = new GeneList(name, Arrays.asList(genes));
        GeneListManager.getInstance().addGeneList(gl);
        session.setCurrentGeneList(gl);

        // Adjust frames
        processFrames(element);
    }

    private void processFrames(Element element) {
        NodeList elements = element.getChildNodes();
        if (elements.getLength() > 0) {
            Map<String, ReferenceFrame> frames = new HashMap();
            for (ReferenceFrame f : FrameManager.getFrames()) {
                frames.put(f.getName(), f);
            }
            List<ReferenceFrame> reorderedFrames = new ArrayList();

            for (int i = 0; i < elements.getLength(); i++) {
                Node childNode = elements.item(i);
                if (childNode.getNodeName().equalsIgnoreCase(SessionElement.FRAME.getText())) {
                    String frameName = getAttribute((Element) childNode, SessionAttribute.NAME.getText());

                    ReferenceFrame f = frames.get(frameName);
                    if (f != null) {
                        reorderedFrames.add(f);
                        try {
                            String chr = getAttribute((Element) childNode, SessionAttribute.CHR.getText());
                            int start = (int) Double.parseDouble(getAttribute((Element) childNode, SessionAttribute.START.getText()));
                            int end = (int) Double.parseDouble(getAttribute((Element) childNode, SessionAttribute.END.getText()));
                            f.setInterval(chr, start, end);
                        } catch (NumberFormatException e) {
                            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                    }

                }
            }
            if (reorderedFrames.size() > 0) {
                FrameManager.setFrames(reorderedFrames);
            }
        }
        IGV.getInstance().resetFrames();
    }

    private void processFilter(Session session, Element element, HashMap additionalInformation) {

        String match = getAttribute(element, SessionAttribute.FILTER_MATCH.getText());
        String showAllTracks = getAttribute(element, SessionAttribute.FILTER_SHOW_ALL_TRACKS.getText());

        String filterName = getAttribute(element, SessionAttribute.NAME.getText());
        TrackFilter filter = new TrackFilter(filterName, null);
        additionalInformation.put(SessionElement.FILTER, filter);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);

        // Save the filter
        session.setFilter(filter);

        // Set filter properties
        if ("all".equalsIgnoreCase(match)) {
            IGV.getInstance().setFilterMatchAll(true);
        } else if ("any".equalsIgnoreCase(match)) {
            IGV.getInstance().setFilterMatchAll(false);
        }

        if ("true".equalsIgnoreCase(showAllTracks)) {
            IGV.getInstance().setFilterShowAllTracks(true);
        } else {
            IGV.getInstance().setFilterShowAllTracks(false);
        }
    }

    private void processFilterElement(Session session, Element element,
                                      HashMap additionalInformation) {

        TrackFilter filter = (TrackFilter) additionalInformation.get(SessionElement.FILTER);
        String item = getAttribute(element, SessionAttribute.ITEM.getText());
        String operator = getAttribute(element, SessionAttribute.OPERATOR.getText());
        String value = getAttribute(element, SessionAttribute.VALUE.getText());
        String booleanOperator = getAttribute(element, SessionAttribute.BOOLEAN_OPERATOR.getText());

        TrackFilterElement trackFilterElement = new TrackFilterElement(filter, item,
                Operator.findEnum(operator), value,
                BooleanOperator.findEnum(booleanOperator));
        filter.add(trackFilterElement);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);
    }

    /**
     * A counter to generate unique panel names.  Needed for backward-compatibility of old session files.
     */
    private int panelCounter = 1;

    private void processPanel(Session session, Element element, HashMap additionalInformation) {
        panelElementPresent = true;
        String panelName = element.getAttribute("name");
        if (panelName == null) {
            panelName = "Panel" + panelCounter++;
        }

        List<Track> panelTracks = new ArrayList();
        NodeList elements = element.getChildNodes();
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            if (childNode.getNodeName().equalsIgnoreCase(SessionElement.DATA_TRACK.getText()) ||  // Is this a track?
                    childNode.getNodeName().equalsIgnoreCase(SessionElement.TRACK.getText())) {

                List<Track> tracks = processTrack(session, (Element) childNode, additionalInformation);
                if (tracks != null) {
                    panelTracks.addAll(tracks);
                }
            } else {
                process(session, childNode, additionalInformation);
            }
        }

        TrackPanel panel = IGV.getInstance().getDataPanel(panelName);
        panel.addTracks(panelTracks);
    }

    private void processPanelLayout(Session session, Element element, HashMap additionalInformation) {

        String nodeName = element.getNodeName();
        String panelName = nodeName;

        NamedNodeMap tNodeMap = element.getAttributes();
        for (int i = 0; i < tNodeMap.getLength(); i++) {
            Node node = tNodeMap.item(i);
            String name = node.getNodeName();
            if (name.equals("dividerFractions")) {
                String value = node.getNodeValue();
                String[] tokens = value.split(",");
                double[] divs = new double[tokens.length];
                try {
                    for (int j = 0; j < tokens.length; j++) {
                        divs[j] = Double.parseDouble(tokens[j]);
                    }
                    session.setDividerFractions(divs);
                }
                catch (NumberFormatException e) {
                    log.error("Error parsing divider locations", e);
                }
            }
        }
    }


    /**
     * Process a track element.  This should return a single track, but could return multiple tracks since the
     * uniqueness of the track id is not enforced.
     *
     * @param session
     * @param element
     * @param additionalInformation
     * @return
     */

    private List<Track> processTrack(Session session, Element element, HashMap additionalInformation) {

        String id = getAttribute(element, SessionAttribute.ID.getText());

        // TODo -- put in utility method, extacts attributes from element **Definitely need to do this
        HashMap<String, String> tAttributes = new HashMap();
        HashMap<String, String> drAttributes = null;

        NamedNodeMap tNodeMap = element.getAttributes();
        for (int i = 0; i < tNodeMap.getLength(); i++) {
            Node node = tNodeMap.item(i);
            String value = node.getNodeValue();
            if (value != null && value.length() > 0) {
                tAttributes.put(node.getNodeName(), value);
            }
        }


        if (element.hasChildNodes()) {
            drAttributes = new HashMap();
            Node childNode = element.getFirstChild();
            Node sibNode = childNode.getNextSibling();
            String sibName = sibNode.getNodeName();
            if (sibName.equals(SessionElement.DATA_RANGE.getText())) {
                NamedNodeMap drNodeMap = sibNode.getAttributes();
                for (int i = 0; i < drNodeMap.getLength(); i++) {
                    Node node = drNodeMap.item(i);
                    String value = node.getNodeValue();
                    if (value != null && value.length() > 0) {
                        drAttributes.put(node.getNodeName(), value);
                    }
                }
            }
        }

        // Get matching tracks.  The trackNameDictionary is used for pre V 2 files, where ID was loosely defined
        List<Track> matchedTracks = trackDictionary.get(id);

        if (matchedTracks == null) {
            log.info("Warning.  No tracks were found with id: " + id + " in session file");
        } else {
            for (final Track track : matchedTracks) {

                // Special case for sequence & gene tracks,  they need to be removed before being placed.
                if (version >= 4 && track == geneTrack || track == seqTrack) {
                    IGV.getInstance().getTrackManager().removeTracks(Arrays.asList(track));
                }

                track.restorePersistentState(tAttributes);
                if (drAttributes != null) {
                    DataRange dr = track.getDataRange();
                    dr.restorePersistentState(drAttributes);
                    track.setDataRange(dr);
                }
            }
            trackDictionary.remove(id);


        }

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);

        return matchedTracks;
    }

    private void processColorScales(Session session, Element element, HashMap additionalInformation) {

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);
    }

    private void processColorScale(Session session, Element element, HashMap additionalInformation) {

        String trackType = getAttribute(element, SessionAttribute.TYPE.getText());
        String value = getAttribute(element, SessionAttribute.VALUE.getText());

        setColorScaleSet(session, trackType, value);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation);
    }

    private void processPreferences(Session session, Element element, HashMap additionalInformation) {

        NodeList elements = element.getChildNodes();
        for (int i = 0; i < elements.getLength(); i++) {
            Node child = elements.item(i);
            if (child.getNodeName().equalsIgnoreCase(SessionElement.PROPERTY.getText())) {
                Element childNode = (Element) child;
                String name = getAttribute(childNode, SessionAttribute.NAME.getText());
                String value = getAttribute(childNode, SessionAttribute.VALUE.getText());
                session.setPreference(name, value);
            }
        }
    }


    /**
     * Process a list of session element nodes.
     *
     * @param session
     * @param elements
     */
    private void process(Session session, NodeList elements, HashMap additionalInformation) {
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            process(session, childNode, additionalInformation);
        }
    }


    /**
     * Reads an xml from an input file and creates DOM document.
     *
     * @param
     * @return
     * @throws ParserConfigurationException
     * @throws IOException
     * @throws SAXException
     */
    private Document createDOMDocumentFromXmlFile(InputStream inputStream)
            throws ParserConfigurationException, IOException, SAXException {
        DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
        Document xmlDocument = documentBuilder.parse(inputStream);
        return xmlDocument;
    }


    public void setColorScaleSet(Session session, String type, String value) {

        if (type == null | value == null) {
            return;
        }

        TrackType trackType = TrackType.OTHER;

        if (TrackType.ALLELE_SPECIFIC_COPY_NUMBER.name().equalsIgnoreCase(type)) {
            trackType = TrackType.ALLELE_SPECIFIC_COPY_NUMBER;
        } else if (TrackType.CHIP.name().equalsIgnoreCase(type)) {
            trackType = TrackType.CHIP;
        } else if (TrackType.COPY_NUMBER.name().equalsIgnoreCase(type)) {
            trackType = TrackType.COPY_NUMBER;
        } else if (TrackType.DNA_METHYLATION.name().equalsIgnoreCase(type)) {
            trackType = TrackType.DNA_METHYLATION;
        } else if (TrackType.OTHER.name().equalsIgnoreCase(type)) {
            trackType = TrackType.OTHER;
        } else if (TrackType.GENE_EXPRESSION.name().equalsIgnoreCase(type)) {
            trackType = TrackType.GENE_EXPRESSION;
        } else if (TrackType.LOH.name().equalsIgnoreCase(type)) {
            trackType = TrackType.LOH;
        } else if (TrackType.MUTATION.name().equalsIgnoreCase(type)) {
            trackType = TrackType.MUTATION;
        } else if (TrackType.PHASTCON.name().equalsIgnoreCase(type)) {
            trackType = TrackType.PHASTCON;
        } else if (TrackType.TILING_ARRAY.name().equalsIgnoreCase(type)) {
            trackType = TrackType.TILING_ARRAY;
        }

        // TODO -- refactor to remove instanceof / cast.  Currently only ContinuousColorScale is handled
        ColorScale colorScale = ColorScaleFactory.getScaleFromString(value);
        if (colorScale instanceof ContinuousColorScale) {
            session.setColorScale(trackType, (ContinuousColorScale) colorScale);
        }

        // ColorScaleFactory.setColorScale(trackType, colorScale);
    }

    private String getAttribute(Element element, String key) {
        String value = element.getAttribute(key);
        if (value != null) {
            if (value.trim().equals("")) {
                value = null;
            }
        }
        return value;
    }


}
