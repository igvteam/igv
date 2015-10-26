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
import org.broad.igv.Globals;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.TrackFilterElement;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.ui.panel.TrackPanelScrollPane;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.FilterElement.BooleanOperator;
import org.broad.igv.util.FilterElement.Operator;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.igv.util.collections.CollUtils;
import org.w3c.dom.*;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ref.WeakReference;
import java.util.*;
import java.util.List;

/**
 * Class to parse an IGV session file
 */
public class IGVSessionReader implements SessionReader {

    private static Logger log = Logger.getLogger(IGVSessionReader.class);
    private static String INPUT_FILE_KEY = "INPUT_FILE_KEY";
    // Temporary values used in processing

    //package-private for unit testing
    Collection<ResourceLocator> dataFiles;
    private Collection<ResourceLocator> missingDataFiles;
    private static Map<String, String> attributeSynonymMap = new HashMap();
    private boolean panelElementPresent = false;
    private int version;

    private IGV igv;

    private static WeakReference<IGVSessionReader> currentReader;


    /**
     * Classes that have been registered for use with JAXB
     */
    private static List<Class> registeredClasses = new ArrayList<Class>();


    /**
     * Map of track id -> track.  It is important to maintain the order in which tracks are added, thus
     * the use of LinkedHashMap. We add tracks here when loaded and remove them when attributes are specified.
     */
    private final Map<String, List<Track>> leftoverTrackDictionary = Collections.synchronizedMap(new LinkedHashMap());

    /**
     * Map of id -> track, for second pass through when tracks reference each other
     */
    private final Map<String, List<Track>> allTracks = Collections.synchronizedMap(new LinkedHashMap<String, List<Track>>());
    private String rootPath;

    public List<Track> getTracksById(String trackId) {
        return allTracks.get(trackId);
    }


    /**
     * Map of full path -> relative path.
     */
    Map<String, String> fullToRelPathMap = new HashMap<String, String>();

    private Track geneTrack = null;
    private Track seqTrack = null;
    private boolean hasTrackElments;

    //Temporary holder for generating tracks
    protected static AbstractTrack nextTrack;

    static {
        attributeSynonymMap.put("DATA FILE", "DATA SET");
        attributeSynonymMap.put("TRACK NAME", "NAME");

        registerClass(AbstractTrack.class);
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
        VISIBLE_ATTRIBUTES("VisibleAttributes"),
        ATTRIBUTE("Attribute"),
        VISIBLE_ATTRIBUTE("VisibleAttribute"),
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
        NEXT_AUTOSCALE_GROUP("nextAutoscaleGroup"),
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
        HAS_GENE_TRACK("hasGeneTrack"),
        HAS_SEQ_TRACK("hasSequenceTrack"),

        //RESOURCE ATTRIBUTES
        PATH("path"),
        INDEX("index"),
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
        ShadeBasesOption shadeBases;
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

    }

    public IGVSessionReader(IGV igv) {
        this.igv = igv;
        currentReader = new WeakReference<IGVSessionReader>(this);
    }


    /**
     * @param inputStream
     * @param session
     * @param sessionPath @return
     * @throws RuntimeException
     */

    public void loadSession(InputStream inputStream, Session session, String sessionPath) {


        log.debug("Load session");


        Document document = null;
        try {
            document = Utilities.createDOMDocumentFromXmlStream(inputStream);
        } catch (Exception e) {
            log.error("Load session error", e);
            throw new RuntimeException(e);
        }

        NodeList tracks = document.getElementsByTagName("Track");
        hasTrackElments = tracks.getLength() > 0;

        HashMap additionalInformation = new HashMap();
        additionalInformation.put(INPUT_FILE_KEY, sessionPath);

        NodeList nodes = document.getElementsByTagName(SessionElement.GLOBAL.getText());
        if (nodes == null || nodes.getLength() == 0) {
            nodes = document.getElementsByTagName(SessionElement.SESSION.getText());
        }

        this.rootPath = sessionPath;
        processRootNode(session, nodes.item(0), additionalInformation, sessionPath);

        // Add tracks not explicitly allocated to panels.  It is legal to define sessions with the Resources
        // section only (no Panel or Track elements).
        addLeftoverTracks(leftoverTrackDictionary.values());

        if (igv != null) {
            if (session.getGroupTracksBy() != null && session.getGroupTracksBy().length() > 0) {
                igv.setGroupByAttribute(session.getGroupTracksBy());
            }

            if (session.isRemoveEmptyPanels()) {
                igv.getMainPanel().removeEmptyDataPanels();
            }

            igv.resetOverlayTracks();
        }

    }


    private void processRootNode(Session session, Node node, HashMap additionalInformation, String rootPath) {

        if ((node == null) || (session == null)) {
            MessageUtils.showMessage("Invalid session file: root node not found");
            return;
        }

        String nodeName = node.getNodeName();
        if (!(nodeName.equalsIgnoreCase(SessionElement.GLOBAL.getText()) || nodeName.equalsIgnoreCase(SessionElement.SESSION.getText()))) {
            MessageUtils.showMessage("Session files must begin with a \"Global\" or \"Session\" element.  Found: " + nodeName);
            return;
        }

        Element element = (Element) node;
        String alternateRootPath = getAttribute(element, SessionAttribute.PATH.getText());

        process(session, node, additionalInformation, rootPath, alternateRootPath);

        String versionString = getAttribute(element, SessionAttribute.VERSION.getText());
        try {
            version = Integer.parseInt(versionString);
        } catch (NumberFormatException e) {
            log.error("Non integer version number in session file: " + versionString);
        }

        // Load the genome, which can be an ID, or a path or URL to a .genome or indexed fasta file.
        String genomeId = getAttribute(element, SessionAttribute.GENOME.getText());
        String hasGeneTrackStr = getAttribute(element, SessionAttribute.HAS_GENE_TRACK.getText());
        boolean hasGeneTrack = true;
        if (hasGeneTrackStr != null) {
            hasGeneTrack = Boolean.parseBoolean(hasGeneTrackStr);
        }
        boolean hasSeqTrack = hasGeneTrack;
        String hasSeqTrackStr = getAttribute(element, SessionAttribute.HAS_SEQ_TRACK.getText());
        if (hasSeqTrackStr != null) {
            hasSeqTrack = Boolean.parseBoolean(hasSeqTrackStr);
        }

        if (genomeId != null && genomeId.length() > 0) {
            if (genomeId.equals(GenomeManager.getInstance().getGenomeId())) {
                // We don't have to reload the genome, but the gene track for the current genome should be restored.
                if (hasGeneTrack || hasSeqTrack) {
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
                    FeatureTrack geneTrack = hasGeneTrack ? genome.getGeneTrack() : null;
                    IGV.getInstance().setGenomeTracks(geneTrack);
                }
            } else {
                // Selecting a genome will actually "reset" the session so we have to
                // save the path and restore it.
                String sessionPath = session.getPath();
                //Loads genome from list, or from server or cache
                igv.selectGenomeFromList(genomeId);
                if (!genomeId.equals(GenomeManager.getInstance().getGenomeId())) {
                    String genomePath = genomeId;
                    if (!ParsingUtils.pathExists(genomePath)) {
                        genomePath = FileUtils.getAbsolutePath(session.getPath(), genomeId);
                    }
                    if (ParsingUtils.pathExists(genomePath)) {
                        try {
                            IGV.getInstance().loadGenome(genomePath, null, hasGeneTrack);
                        } catch (IOException e) {
                            throw new RuntimeException("Error loading genome: " + genomeId);
                        }
                    } else {
                        MessageUtils.showMessage("Warning: Could not locate genome: " + genomeId);
                    }
                }
                session.setPath(sessionPath);
            }
        }

        //TODO Remove these nearly identical if/then statements

        if (!hasGeneTrack && igv.hasGeneTrack()) {
            //Need to remove gene track if it was loaded because it's not supposed to be in the session
            igv.removeTracks(Arrays.<Track>asList(GenomeManager.getInstance().getCurrentGenome().getGeneTrack()));
            geneTrack = null;
        } else {
            //For later lookup and to prevent dual adding, we keep a reference to the gene track
            geneTrack = GenomeManager.getInstance().getCurrentGenome().getGeneTrack();
            if (geneTrack != null) {
                allTracks.put(geneTrack.getId(), Arrays.asList(geneTrack));
            }
        }

        SequenceTrack tmpSeqTrack = igv.getSequenceTrack();
        if (hasSeqTrack && !igv.hasSequenceTrack()) {
            //This will create a sequence track
            IGV.getInstance().setGenomeTracks(null);
        } else if (!hasSeqTrack && igv.hasSequenceTrack()) {
            //Need to remove seq track if it was loaded because it's not supposed to be in the session
            igv.removeTracks(Arrays.<Track>asList(tmpSeqTrack));
            seqTrack = null;
        } else {
            //For later lookup and to prevent dual adding, we keep a reference to the sequence track
            seqTrack = tmpSeqTrack;
            if (seqTrack != null) {
                allTracks.put(seqTrack.getId(), Arrays.asList(seqTrack));
            }
        }

        session.setLocus(getAttribute(element, SessionAttribute.LOCUS.getText()));
        session.setGroupTracksBy(getAttribute(element, SessionAttribute.GROUP_TRACKS_BY.getText()));

        String nextAutoscaleGroup = getAttribute(element, SessionAttribute.NEXT_AUTOSCALE_GROUP.getText());
        if(nextAutoscaleGroup != null) {
            try {
                session.setNextAutoscaleGroup(Integer.parseInt(nextAutoscaleGroup));
            } catch (NumberFormatException e) {
                log.error("Error setting next autoscale group", e);
            }
        }

        String removeEmptyTracks = getAttribute(element, "removeEmptyTracks");
        if (removeEmptyTracks != null) {
            try {
                Boolean b = Boolean.parseBoolean(removeEmptyTracks);
                session.setRemoveEmptyPanels(b);
            } catch (Exception e) {
                log.error("Error parsing removeEmptyTracks string: " + removeEmptyTracks, e);
            }
        }

        session.setVersion(version);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);

        // ReferenceFrame.getInstance().invalidateLocationScale();
    }

    //TODO Check to make sure tracks are not being created twice
    //TODO -- DONT DO THIS FOR NEW SESSIONS
    private void addLeftoverTracks(Collection<List<Track>> tmp) {
        Map<String, TrackPanel> trackPanelCache = new HashMap();
        if (version < 3 || !panelElementPresent) {
            log.debug("Adding \"leftover\" tracks");

            //For resetting track panels later
            List<Map<TrackPanelScrollPane, Integer>> trackPanelAttrs = null;
            if (IGV.hasInstance()) {
                trackPanelAttrs = IGV.getInstance().getTrackPanelAttrs();
            }

            for (List<Track> tracks : tmp) {
                for (Track track : tracks) {
                    if (track != geneTrack && track != seqTrack && track.getResourceLocator() != null) {

                        TrackPanel panel = trackPanelCache.get(track.getResourceLocator().getPath());
                        if (panel == null) {
                            panel = IGV.getInstance().getPanelFor(track.getResourceLocator());
                            trackPanelCache.put(track.getResourceLocator().getPath(), panel);
                        }
                        panel.addTrack(track);
                    }
                }
            }

            if (IGV.hasInstance()) {
                IGV.getInstance().resetPanelHeights(trackPanelAttrs.get(0), trackPanelAttrs.get(1));
            }
        }

    }


    /**
     * Process a single session element node.
     *
     * @param session
     * @param element
     */
    private void process(Session session, Node element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        if ((element == null) || (session == null)) {
            return;
        }

        String nodeName = element.getNodeName();

        if (nodeName.equalsIgnoreCase(SessionElement.RESOURCES.getText()) ||
                nodeName.equalsIgnoreCase(SessionElement.FILES.getText())) {
            processResources(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.RESOURCE.getText()) ||
                nodeName.equalsIgnoreCase(SessionElement.DATA_FILE.getText())) {
            processResource(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.REGIONS.getText())) {
            processRegions(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.REGION.getText())) {
            processRegion(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.GENE_LIST.getText())) {
            processGeneList(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.FILTER.getText())) {
            processFilter(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.FILTER_ELEMENT.getText())) {
            processFilterElement(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.COLOR_SCALES.getText())) {
            processColorScales(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.COLOR_SCALE.getText())) {
            processColorScale(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.PREFERENCES.getText())) {
            processPreferences(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.DATA_TRACKS.getText()) ||
                nodeName.equalsIgnoreCase(SessionElement.FEATURE_TRACKS.getText()) ||
                nodeName.equalsIgnoreCase(SessionElement.PANEL.getText())) {
            processPanel(session, (Element) element, additionalInformation, rootPath, alternateRootPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.PANEL_LAYOUT.getText())) {
            processPanelLayout(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.HIDDEN_ATTRIBUTES.getText())) {
            processHiddenAttributes(session, (Element) element, additionalInformation);
        } else if (nodeName.equalsIgnoreCase(SessionElement.VISIBLE_ATTRIBUTES.getText())) {
            processVisibleAttributes(session, (Element) element, additionalInformation);
        }

    }

    private void processResources(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {
        dataFiles = new ArrayList();
        missingDataFiles = new ArrayList();
        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);

        if (missingDataFiles.size() > 0) {
            StringBuffer message = new StringBuffer();
            message.append("<html>The following data file(s) could not be located.<ul>");
            for (ResourceLocator file : missingDataFiles) {
                if (file.getDBUrl() == null) {
                    message.append("<li>");
                    message.append(file.getPath());
                    message.append("</li>");
                } else {
                    message.append("<li>Server: ");
                    message.append(file.getDBUrl());
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

            final List<String> errors = new ArrayList<String>();

            // Load files concurrently -- TODO, put a limit on # of threads?
            List<Thread> threads = new ArrayList(dataFiles.size());
            long t0 = System.currentTimeMillis();
            int i = 0;
            List<Runnable> synchronousLoads = new ArrayList<Runnable>();
            for (final ResourceLocator locator : dataFiles) {

                final String suppliedPath = locator.getPath();
                final String relPath = fullToRelPathMap.get(suppliedPath);

                Runnable runnable = new Runnable() {
                    public void run() {
                        List<Track> tracks = null;
                        try {
                            tracks = igv.load(locator);
                            for (Track track : tracks) {
                                if (track == null) {
                                    log.info("Null track for resource " + locator.getPath());
                                    continue;
                                }

                                String id = track.getId();
                                if (id == null) {
                                    log.info("Null track id for resource " + locator.getPath());
                                    continue;
                                }

                                // if (relPath != null) {
                                //     id = id.replace(suppliedPath, relPath);
                                // }

                                List<Track> trackList = leftoverTrackDictionary.get(id);
                                if (trackList == null) {
                                    trackList = new ArrayList();
                                    leftoverTrackDictionary.put(id, trackList);
                                    allTracks.put(id, trackList);
                                }
                                trackList.add(track);
                            }
                        } catch (Exception e) {
                            log.error("Error loading resource " + locator.getPath(), e);
                            String ms = "<b>" + locator.getPath() + "</b><br>&nbs;p&nbsp;" + e.toString() + "<br>";
                            errors.add(ms);
                        }
                    }
                };

                boolean isAlignment = locator.getPath().endsWith(".bam") || locator.getPath().endsWith(".entries") ||
                        locator.getPath().endsWith(".sam");


                // Run synchronously if in batch mode or if there are no "track" elments, or if this is an alignment file
                // EVERYTHING IS RUN SYNCHRONOUSLY FOR NOW UNTIL WE CAN FIGURE OUT WHAT TO DO TO PREVENT MULTIPLE
                // AUTHENTICATION DIALOGS
                if (isAlignment || Globals.isBatch() || !hasTrackElments) {
                    synchronousLoads.add(runnable);
                } else {
                    Thread t = new Thread(runnable);
                    threads.add(t);
                    t.start();
                }
                i++;
            }
            // Wait for all threads to complete
            for (Thread t : threads) {
                try {
                    t.join();
                } catch (InterruptedException ignore) {
                }
            }


            // Now load data that must be loaded synchronously
            for (Runnable runnable : synchronousLoads) {
                runnable.run();
            }

            long dt = System.currentTimeMillis() - t0;
            log.debug("Total load time = " + dt);

            if (errors.size() > 0) {
                StringBuffer buf = new StringBuffer();
                buf.append("<html>Errors were encountered loading the session:<br>");
                for (String msg : errors) {
                    buf.append(msg);
                }
                MessageUtils.showMessage(buf.toString());
            }

        }
        dataFiles = null;
    }

    /**
     * Load a single resource.
     * <p/>
     * Package private for unit testing
     *
     * @param session
     * @param element
     * @param additionalInformation
     */
    void processResource(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        String nodeName = element.getNodeName();
        boolean oldSession = nodeName.equals(SessionElement.DATA_FILE.getText());

        String label = getAttribute(element, SessionAttribute.LABEL.getText());
        String name = getAttribute(element, SessionAttribute.NAME.getText());
        String sampleId = getAttribute(element, SessionAttribute.SAMPLE_ID.getText());
        String description = getAttribute(element, SessionAttribute.DESCRIPTION.getText());
        String type = getAttribute(element, SessionAttribute.TYPE.getText());
        String coverage = getAttribute(element, SessionAttribute.COVERAGE.getText());
        String trackLine = getAttribute(element, SessionAttribute.TRACK_LINE.getText());
        String colorString = getAttribute(element, SessionAttribute.COLOR.getText());
        String index = getAttribute(element, SessionAttribute.INDEX.getText());

        // URL to a database or webservice
        String serverURL = getAttribute(element, SessionAttribute.SERVER_URL.getText());

        // Older sessions used the "name" attribute for the path.
        String path = getAttribute(element, SessionAttribute.PATH.getText());

        if (oldSession && name != null) {
            path = name;
            int idx = name.lastIndexOf("/");
            if (idx > 0 && idx + 1 < name.length()) {
                name = name.substring(idx + 1);
            }
        }

        if (rootPath == null) {
            log.error("Null root path -- this is not expected");
            MessageUtils.showMessage("Unexpected error loading session: null root path");
            return;
        }

        String absolutePath = getAbsolutePath(path, rootPath, alternateRootPath);

        fullToRelPathMap.put(absolutePath, path);

        ResourceLocator resourceLocator = new ResourceLocator(serverURL, absolutePath);

        if (index != null) resourceLocator.setIndexPath(index);

        if (coverage != null) {
            String absoluteCoveragePath = getAbsolutePath(coverage, rootPath, alternateRootPath);
            resourceLocator.setCoverage(absoluteCoveragePath);
        }


        String url = getAttribute(element, SessionAttribute.URL.getText());
        if (url == null) {
            url = getAttribute(element, SessionAttribute.FEATURE_URL.getText());
        }
        resourceLocator.setFeatureInfoURL(url);

        String infolink = getAttribute(element, SessionAttribute.HYPERLINK.getText());
        if (infolink == null) {
            infolink = getAttribute(element, SessionAttribute.INFOLINK.getText());
        }
        resourceLocator.setTrackInforURL(infolink);


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

        dataFiles.add(resourceLocator);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);

    }

    private String getAbsolutePath(String path, String rootPath, String alternateRootPath) {

        String absolutePath = FileUtils.getAbsolutePath(path, rootPath);

        // Session file might have been moved, try to correct if the path is a local file (not an http or ftp url)
        if (alternateRootPath != null && !FileUtils.isRemote(absolutePath) && !(new File(absolutePath).exists())) {
            absolutePath = FileUtils.getAbsolutePath(path, alternateRootPath);
        }

        return absolutePath;
    }

    private void processRegions(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        session.clearRegionsOfInterest();
        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);
    }

    private void processRegion(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        String chromosome = getAttribute(element, SessionAttribute.CHROMOSOME.getText());
        String start = getAttribute(element, SessionAttribute.START_INDEX.getText());
        String end = getAttribute(element, SessionAttribute.END_INDEX.getText());
        String description = getAttribute(element, SessionAttribute.DESCRIPTION.getText());

        RegionOfInterest region = new RegionOfInterest(chromosome, new Integer(start), new Integer(end), description);
        IGV.getInstance().addRegionOfInterest(region);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);
    }


    private void processHiddenAttributes(Session session, Element element, HashMap additionalInformation) {

//        session.clearRegionsOfInterest();
        NodeList elements = element.getChildNodes();
        if (elements.getLength() > 0) {
            Set<String> attributes = new HashSet();
            for (int i = 0; i < elements.getLength(); i++) {
                Node childNode = elements.item(i);
                if (childNode.getNodeName().equals(IGVSessionReader.SessionElement.ATTRIBUTE.getText())) {
                    attributes.add(((Element) childNode).getAttribute(IGVSessionReader.SessionAttribute.NAME.getText()));
                }
            }
            session.setHiddenAttributes(attributes);
        }
    }


    /**
     * For backward compatibility
     *
     * @param session
     * @param element
     * @param additionalInformation
     */
    private void processVisibleAttributes(Session session, Element element, HashMap additionalInformation) {

//        session.clearRegionsOfInterest();
        NodeList elements = element.getChildNodes();
        if (elements.getLength() > 0) {
            Set<String> visibleAttributes = new HashSet();
            for (int i = 0; i < elements.getLength(); i++) {
                Node childNode = elements.item(i);
                if (childNode.getNodeName().equals(IGVSessionReader.SessionElement.VISIBLE_ATTRIBUTE.getText())) {
                    visibleAttributes.add(((Element) childNode).getAttribute(IGVSessionReader.SessionAttribute.NAME.getText()));
                }
            }

            final List<String> attributeNames = AttributeManager.getInstance().getAttributeNames();
            Set<String> hiddenAttributes = new HashSet<String>(attributeNames);
            hiddenAttributes.removeAll(visibleAttributes);
            session.setHiddenAttributes(hiddenAttributes);

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
                            final String startString =
                                    getAttribute((Element) childNode, SessionAttribute.START.getText()).replace(",", "");
                            final String endString =
                                    getAttribute((Element) childNode, SessionAttribute.END.getText()).replace(",", "");
                            int start = ParsingUtils.parseInt(startString);
                            int end = ParsingUtils.parseInt(endString);
                            org.broad.igv.feature.Locus locus = new Locus(chr, start, end);
                            f.jumpTo(locus);
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

    private void processFilter(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        String match = getAttribute(element, SessionAttribute.FILTER_MATCH.getText());
        String showAllTracks = getAttribute(element, SessionAttribute.FILTER_SHOW_ALL_TRACKS.getText());

        String filterName = getAttribute(element, SessionAttribute.NAME.getText());
        TrackFilter filter = new TrackFilter(filterName, null);
        additionalInformation.put(SessionElement.FILTER, filter);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);

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
                                      HashMap additionalInformation, String rootPath, String alternateRootPath) {

        TrackFilter filter = (TrackFilter) additionalInformation.get(SessionElement.FILTER);
        String item = getAttribute(element, SessionAttribute.ITEM.getText());
        String operator = getAttribute(element, SessionAttribute.OPERATOR.getText());
        String value = getAttribute(element, SessionAttribute.VALUE.getText());
        String booleanOperator = getAttribute(element, SessionAttribute.BOOLEAN_OPERATOR.getText());

        Operator opEnum = CollUtils.findValueOf(Operator.class, operator);
        BooleanOperator boolEnum = BooleanOperator.valueOf(booleanOperator.toUpperCase());
        TrackFilterElement trackFilterElement = new TrackFilterElement(filter, item,
                opEnum, value, boolEnum);
        filter.add(trackFilterElement);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);
    }

    /**
     * A counter to generate unique panel names.  Needed for backward-compatibility of old session files.
     */
    private int panelCounter = 1;

    /**
     * @param node
     * @return Whether this node is a track
     */
    private static boolean nodeIsTrack(Node node) {
        return node.getNodeName() != null &&
                (node.getNodeName().equalsIgnoreCase(SessionElement.DATA_TRACK.getText()) ||
                        node.getNodeName().equalsIgnoreCase(SessionElement.TRACK.getText()));
    }

    private void processPanel(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {
        panelElementPresent = true;
        String panelName = element.getAttribute("name");
        if (panelName == null) {
            panelName = "Panel" + panelCounter++;
        }

        List<Track> panelTracks = new ArrayList();
        NodeList elements = element.getChildNodes();
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            if (nodeIsTrack(childNode)) {
                List<Track> tracks = processTrack(session, (Element) childNode, additionalInformation, rootPath, alternateRootPath);
                if (tracks != null) {
                    panelTracks.addAll(tracks);
                }
            } else {
                process(session, childNode, additionalInformation, rootPath, alternateRootPath);
            }
        }

        //We make a second pass through, resolving references to tracks which may have been processed afterwards.
        //For instance if Track 2 referenced Track 4
        //TODO Make this less hacky
        for (Track track : panelTracks) {
            if (track instanceof FeatureTrack) {
                FeatureTrack featureTrack = (FeatureTrack) track;
                featureTrack.updateTrackReferences(panelTracks);
            } else if (track instanceof DataSourceTrack) {
                DataSourceTrack dataTrack = (DataSourceTrack) track;
                dataTrack.updateTrackReferences(panelTracks);
            }
        }

        TrackPanel panel = IGV.getInstance().getTrackPanel(panelName);
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
                } catch (NumberFormatException e) {
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

    private List<Track> processTrack(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        String id = getAttribute(element, SessionAttribute.ID.getText());

        // Get matching tracks.
        List<Track> matchedTracks = allTracks.get(id);

        if (matchedTracks == null) {
            //Try creating an "absolute" path for the id
            matchedTracks = allTracks.get(getAbsolutePath(id, rootPath, alternateRootPath));
        }

        if (matchedTracks == null) {

            log.info("Warning.  No tracks were found with id: " + id + " in session file");
            String className = getAttribute(element, "clazz");

            //We try anyway, some tracks can be reconstructed without a resource element
            //They must have children in that case though, either a source (analysis tracks)
            //or another track (MergedTracks)
            try {
                if (className != null && element.hasChildNodes()) {
                    Track track = null;
                    if (className.contains("FeatureTrack") || className.contains("DataSourceTrack")) {
                        Class clazz = Class.forName(className);
                        Unmarshaller u = getJAXBContext().createUnmarshaller();
                        track = unmarshalTrackElement(u, element, null, clazz);
                        matchedTracks = new ArrayList<Track>(Arrays.asList(track));
                    } else if (className.contains("MergedTracks")) {
                        List<Track> childTracks = processChildTracks(session, element,
                                additionalInformation, rootPath, alternateRootPath);
                        List<DataTrack> memberTracks = new ArrayList<DataTrack>(childTracks.size());

                        for (Track aTrack : childTracks) {
                            memberTracks.add((DataTrack) aTrack);
                        }
                        track = new MergedTracks(id,
                                getAttribute(element, SessionAttribute.NAME.getText()), memberTracks);
                        matchedTracks = Arrays.asList(track);
                    }
                    if (track != null) {
                        allTracks.put(track.getId(), matchedTracks);
                    }
                }
            } catch (JAXBException e) {
                //pass
            } catch (ClassNotFoundException e) {
                //pass
            }


        } else {

            try {
                Unmarshaller u = getJAXBContext().createUnmarshaller();
                for (final Track track : matchedTracks) {

                    // Special case for sequence & gene tracks, they need to be removed before being placed.
                    if (igv != null && version >= 4 && (track == geneTrack || track == seqTrack)) {
                        igv.removeTracks(Arrays.asList(track));
                    }
                    unmarshalTrackElement(u, element, (AbstractTrack) track);
                }

            } catch (JAXBException e) {
                throw new RuntimeException(e);
            }
            leftoverTrackDictionary.remove(id);
        }

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);

        return matchedTracks;
    }

    /**
     * Recursively loop through children of {@code element}, and process them as tracks iff they are determined
     * to be so
     *
     * @param session
     * @param element
     * @param additionalInformation
     * @param rootPath
     * @return List of processed tracks.
     */
    private List<Track> processChildTracks(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        NodeList memberTrackNodes = element.getChildNodes();
        List<Track> memberTracks = new ArrayList<Track>(memberTrackNodes.getLength());
        for (int index = 0; index < memberTrackNodes.getLength(); index++) {
            Node memberNode = memberTrackNodes.item(index);
            if (nodeIsTrack(memberNode)) {
                List<Track> addedTracks = processTrack(session, (Element) memberNode, additionalInformation, rootPath, alternateRootPath);
                if (addedTracks != null) {
                    memberTracks.addAll(addedTracks);
                }

            }
        }
        return memberTracks;
    }

    private static void setNextTrack(AbstractTrack track) {
        nextTrack = track;
    }

    /**
     * Used for unmarshalling track; JAXB needs a static no-arg factory method
     *
     * @return
     */
    public static AbstractTrack getNextTrack() {
        return nextTrack;
    }

    /**
     * Unmarshal element into specified class
     *
     * @param u
     * @param e
     * @param track
     * @return
     * @throws JAXBException
     */
    protected Track unmarshalTrackElement(Unmarshaller u, Element e, AbstractTrack track) throws JAXBException {
        return unmarshalTrackElement(u, e, track, track.getClass());
    }

    /**
     * @param u
     * @param element
     * @param track      The track into which to unmarshal. Can be null if the relevant static factory method can handle
     *                   creating a new instance
     * @param trackClass Class of track to use for unmarshalling
     * @return The unmarshalled track
     * @throws JAXBException
     */
    protected Track unmarshalTrackElement(Unmarshaller u, Element element, AbstractTrack track, Class trackClass) throws JAXBException {
        AbstractTrack ut;

        synchronized (IGVSessionReader.class) {
            setNextTrack(track);
            ut = unmarshalTrack(u, element, trackClass, trackClass);
        }
        ut.restorePersistentState(element, version);
        return ut;
    }

    private void processColorScales(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);
    }

    private void processColorScale(Session session, Element element, HashMap additionalInformation, String rootPath, String alternateRootPath) {

        String trackType = getAttribute(element, SessionAttribute.TYPE.getText());
        String value = getAttribute(element, SessionAttribute.VALUE.getText());

        setColorScaleSet(session, trackType, value);

        NodeList elements = element.getChildNodes();
        process(session, elements, additionalInformation, rootPath, alternateRootPath);
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
    private void process(Session session, NodeList elements, HashMap additionalInformation, String rootPath, String alternateRootPath) {
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            process(session, childNode, additionalInformation, rootPath, alternateRootPath);
        }
    }

    public void setColorScaleSet(Session session, String type, String value) {

        if (type == null | value == null) {
            return;
        }

        TrackType trackType = CollUtils.valueOf(TrackType.class, type.toUpperCase(), TrackType.OTHER);

        // TODO -- refactor to remove instanceof / cast.  Currently only ContinuousColorScale is handled
        ColorScale colorScale = ColorScaleFactory.getScaleFromString(value);
        if (colorScale instanceof ContinuousColorScale) {
            session.setColorScale(trackType, (ContinuousColorScale) colorScale);
        }

        // ColorScaleFactory.setColorScale(trackType, colorScale);
    }

    private static String getAttribute(Element element, String key) {
        String value = element.getAttribute(key);
        if (value != null) {
            if (value.trim().equals("")) {
                value = null;
            }
        }
        return value;
    }


    private static JAXBContext jc = null;

    public static synchronized JAXBContext getJAXBContext() throws JAXBException {
        if (jc == null) {
            jc = JAXBContext.newInstance(registeredClasses.toArray(new Class[0]), new HashMap<String, Object>());
        }
        return jc;
    }

    /**
     * Register this class with JAXB, so it can be saved and restored to a session.
     * The class must conform the JAXBs requirements (e.g. no-arg constructor or factory method)
     *
     * @param clazz
     */
    //@api
    public static synchronized void registerClass(Class clazz) {
        registeredClasses.add(clazz);
        jc = null;
    }


    /**
     * Unmarshal node. We first attempt to unmarshal into the specified {@code clazz}
     * if that fails, we try the superclass, and so on up.
     *
     * @param node
     * @param unmarshalClass Class to which to use for unmarshalling
     * @param firstClass     The first class used for invocation. For helpful error message only
     * @return
     */
    public static AbstractTrack unmarshalTrack(Unmarshaller u, Node node, Class unmarshalClass, Class firstClass) throws JAXBException {

        if (unmarshalClass == null || unmarshalClass.equals(Object.class)) {
            throw new JAXBException(firstClass + " and none of its superclasses are known");
        }

        if (AbstractTrack.knownUnknownTrackClasses.contains(unmarshalClass)) {
            return unmarshalTrack(u, node, firstClass, unmarshalClass.getSuperclass());
        }

        JAXBElement el;
        try {
            el = u.unmarshal(node, unmarshalClass);
        } catch (JAXBException e) {
            AbstractTrack.knownUnknownTrackClasses.add(unmarshalClass);
            return unmarshalTrack(u, node, firstClass, unmarshalClass.getSuperclass());
        }
        return (AbstractTrack) el.getValue();
    }

    /**
     * Uses {@link #currentReader} to lookup matching tracks by id, or
     * searches allTracks if sessionReader is null
     *
     * @param trackId
     * @param allTracks
     * @return
     */
    public static Track getMatchingTrack(String trackId, List<Track> allTracks) {
        IGVSessionReader reader = currentReader.get();
        List<Track> matchingTracks;
        if (reader != null) {
            matchingTracks = reader.getTracksById(trackId);
        } else {
            if (allTracks == null)
                throw new IllegalStateException("No session reader and no tracks to search to resolve Track references");
            matchingTracks = new ArrayList<Track>();
            for (Track track : allTracks) {
                if (trackId.equals(track.getId())) {
                    matchingTracks.add(track);
                    break;
                }
            }
        }
        if (matchingTracks == null || matchingTracks.size() == 0) {
            //Either the session file is corrupted, or we just haven't loaded the relevant track yet
            return null;
        } else if (matchingTracks.size() >= 2) {
            log.debug("Found multiple tracks with id  " + trackId + ", using the first");
        }
        return matchingTracks.get(0);
    }

}
