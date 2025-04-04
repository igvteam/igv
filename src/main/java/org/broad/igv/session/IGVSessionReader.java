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

import org.broad.igv.Globals;
import org.broad.igv.bedpe.InteractionTrack;
import org.broad.igv.data.CombinedDataSource;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.basepair.BasePairTrack;
import org.broad.igv.feature.dsi.DSITrack;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.sprite.ClusterTrack;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.maf.MultipleAlignmentTrack;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.sam.EWigTrack;
import org.broad.igv.sam.SpliceJunctionTrack;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.TrackFilterElement;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.ui.panel.TrackPanelScrollPane;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.igv.util.FilterElement.BooleanOperator;
import org.broad.igv.util.FilterElement.Operator;
import org.broad.igv.util.collections.CollUtils;
import org.broad.igv.variant.VariantTrack;
import org.w3c.dom.*;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ref.WeakReference;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Class to parse an IGV session file
 */
public class IGVSessionReader implements SessionReader {

    private static Logger log = LogManager.getLogger(IGVSessionReader.class);
    private static Map<String, String> attributeSynonymMap = new HashMap();
    private static WeakReference<IGVSessionReader> currentReader;

    private IGV igv;
    private int version;
    private String genomePath;
    private Collection<ResourceLocator> dataFiles;
    private Collection<ResourceLocator> missingDataFiles;
    private boolean panelElementPresent = false;    // Flag indicating if "Panel" sections are present
    private TrackFilter filter;  // There is a single TrackFilter object per session, usually null

    /**
     * List of combined data source tracks.  Processing of combined data sources has to be deferred until all tracks
     * are loaded
     */
    private final List<Pair<CombinedDataTrack, Element>> combinedDataSourceTracks = new ArrayList<>();

    private Set<Track> allocatedToPanel = new LinkedHashSet<>();;  // List of tracks allocated to panels, if Panel elements are present.

    /**
     * Map of id -> track, for second pass through when tracks reference each other
     */
    private final Map<String, List<Track>> allTracks = Collections.synchronizedMap(new LinkedHashMap<>());

    private final Set<String> erroredResources = Collections.synchronizedSet(new HashSet<>());

    private boolean hasTrackElments;

    static {
        attributeSynonymMap.put("DATA FILE", "DATA SET");
        attributeSynonymMap.put("TRACK NAME", "NAME");
    }


    public IGVSessionReader(IGV igv) {
        this.igv = igv;
        currentReader = new WeakReference<IGVSessionReader>(this);
    }


    /**
     * @param inputStream
     * @param session
     * @param sessionPath -- path to the session file.  This can be null
     * @throws RuntimeException
     */

    public void loadSession(InputStream inputStream, Session session, String sessionPath) {

        Document document;
        try {
            document = Utilities.createDOMDocumentFromXmlStream(inputStream);
        } catch (Exception e) {
            log.error("Load session error", e);
            throw new RuntimeException(e);
        }

        NodeList trackElements = document.getElementsByTagName("Track");
        hasTrackElments = trackElements.getLength() > 0;

        //  Find the root node, either "session" or "global
        NodeList nodes = document.getElementsByTagName(SessionElement.SESSION);
        if (nodes == null || nodes.getLength() == 0) {
            nodes = document.getElementsByTagName(SessionElement.GLOBAL);
        }
        Node rootNode = nodes.item(0);

        // Walk tree processing all nodes, starting with the root
        processRootNode(session, rootNode, sessionPath);

        processCombinedDataSourceTracks();

        // Add tracks not explicitly allocated to panels.   This can happen if a track resource path changes after
        // session created, and session is hand-editted.  It can also happen if the annotation paths for a genome change
        // after session creation.

        List<Track> unallocatedTracks = new ArrayList<>();
        for (List<Track> tracks : allTracks.values()) {
            for (Track t : tracks) {
                if (allocatedToPanel == null ||!allocatedToPanel.contains(t)) {
                    unallocatedTracks.add(t);
                }
            }
        }
        addUnallocatedTracks(unallocatedTracks);


        if (session.getGroupTracksBy() != null && session.getGroupTracksBy().length() > 0) {
            igv.setGroupByAttribute(session.getGroupTracksBy());
        }

        if (session.isRemoveEmptyPanels()) {
            igv.getMainPanel().removeEmptyDataPanels();
        }

        igv.resetOverlayTracks();

    }

    /**
     * The main entry point
     *
     * @param session
     * @param node
     * @param sessionPath
     */

    private void processRootNode(Session session, Node node, String sessionPath) {
        if (node == null || session == null) {
            MessageUtils.showMessage("Invalid session file: root node not found");
            return;
        }

        String nodeName = node.getNodeName();
        if (!nodeName.equalsIgnoreCase(SessionElement.GLOBAL) && !nodeName.equalsIgnoreCase(SessionElement.SESSION)) {
            MessageUtils.showMessage("Session files must begin with a \"Global\" or \"Session\" element. Found: " + nodeName);
            return;
        }

        Element rootElement = (Element) node;
        String versionString = getAttribute(rootElement, SessionAttribute.VERSION);
        try {
            version = Integer.parseInt(versionString);
        } catch (NumberFormatException e) {
            log.error("Non integer version number in session file: " + versionString);
        }

        String genomeId = getAttribute(rootElement, SessionAttribute.GENOME);
        if (genomeId != null && !genomeId.isEmpty()) {
            handleGenomeLoading(sessionPath, genomeId);
        }

        session.setLocus(getAttribute(rootElement, SessionAttribute.LOCUS));
        session.setGroupTracksBy(getAttribute(rootElement, SessionAttribute.GROUP_TRACKS_BY));
        setNextAutoscaleGroup(session, getAttribute(rootElement, SessionAttribute.NEXT_AUTOSCALE_GROUP));
        setRemoveEmptyTracks(session, getAttribute(rootElement, "removeEmptyTracks"));

        session.setVersion(version);
        process(session, rootElement.getChildNodes(), sessionPath);
    }

    private void handleGenomeLoading(String sessionPath, String genomeId) {
        if (genomeId.equals(GenomeManager.getInstance().getGenomeId())) {
            igv.resetSession(sessionPath);
            GenomeManager.getInstance().restoreGenomeTracks(GenomeManager.getInstance().getCurrentGenome());
        } else {
            try {
                GenomeListItem item = GenomeListManager.getInstance().getGenomeListItem(genomeId);
                genomePath = (item != null) ? item.getPath() : getAbsolutePath(genomeId, sessionPath);
                GenomeManager.getInstance().loadGenome(genomePath);
            } catch (IOException e) {
                MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                log.error("Error loading genome: " + genomeId, e);
            }
        }
    }

    private void setNextAutoscaleGroup(Session session, String nextAutoscaleGroup) {
        if (nextAutoscaleGroup != null) {
            try {
                session.setNextAutoscaleGroup(Integer.parseInt(nextAutoscaleGroup));
            } catch (NumberFormatException e) {
                log.error("Error setting next autoscale group", e);
            }
        }
    }

    private void setRemoveEmptyTracks(Session session, String removeEmptyTracks) {
        if (removeEmptyTracks != null) {
            try {
                session.setRemoveEmptyPanels(Boolean.parseBoolean(removeEmptyTracks));
            } catch (Exception e) {
                log.error("Error parsing removeEmptyTracks string: " + removeEmptyTracks, e);
            }
        }
    }

    private void addUnallocatedTracks(List<Track> tracks) {
        Map<String, TrackPanel> trackPanelCache = new HashMap();

        log.debug("Adding \"leftover\" tracks");
        //For resetting track panels later
        List<Map<TrackPanelScrollPane, Integer>> trackPanelAttrs = null;
        trackPanelAttrs = igv.getTrackPanelAttrs();

        for (Track track : tracks) {
            if (track.getResourceLocator() != null) {
                TrackPanel panel = trackPanelCache.get(track.getResourceLocator().getPath());
                if (panel == null) {
                    panel = igv.getPanelFor(track);
                    trackPanelCache.put(track.getResourceLocator().getPath(), panel);
                }
                panel.addTrack(track);
            }
        }

        igv.resetPanelHeights(trackPanelAttrs.get(0), trackPanelAttrs.get(1));
    }


    /**
     * Process a single session element node.
     *
     * @param session
     * @param element
     */
    private void process(Session session, Node element, String sessionPath) {

        if ((element == null) || (session == null)) {
            return;
        }

        String nodeName = element.getNodeName();
        if (nodeName.equalsIgnoreCase(SessionElement.RESOURCES) ||
                nodeName.equalsIgnoreCase(SessionElement.FILES)) {
            processResources(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.RESOURCE) ||
                nodeName.equalsIgnoreCase(SessionElement.DATA_FILE)) {
            processResource(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.REGIONS)) {
            processRegions(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.REGION)) {
            processRegion(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.GENE_LIST)) {
            processGeneList(session, (Element) element);
        } else if (nodeName.equalsIgnoreCase(SessionElement.FILTER)) {
            processFilter(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.FILTER_ELEMENT)) {
            processFilterElement(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.COLOR_SCALES)) {
            processColorScales(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.COLOR_SCALE)) {
            processColorScale(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.PREFERENCES)) {
            processPreferences(session, (Element) element);
        } else if (nodeName.equalsIgnoreCase(SessionElement.DATA_TRACKS) ||
                nodeName.equalsIgnoreCase(SessionElement.FEATURE_TRACKS) ||
                nodeName.equalsIgnoreCase(SessionElement.PANEL)) {
            processPanel(session, (Element) element, sessionPath);
        } else if (nodeName.equalsIgnoreCase(SessionElement.PANEL_LAYOUT)) {
            processPanelLayout(session, (Element) element);
        } else if (nodeName.equalsIgnoreCase(SessionElement.HIDDEN_ATTRIBUTES)) {
            processHiddenAttributes(session, (Element) element);
        } else if (nodeName.equalsIgnoreCase(SessionElement.VISIBLE_ATTRIBUTES)) {
            processVisibleAttributes(session, (Element) element);
        }

    }

    private void processResources(Session session, Element element, String sessionPath) {

        dataFiles = new ArrayList();
        missingDataFiles = new ArrayList();
        NodeList elements = element.getChildNodes();

        process(session, elements, sessionPath);

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

        // Filter data files that are included in genome annotations.  This is to work around a bug introduced with
        // json genomes where the genome annotation resources were (doubly) included as session
        List<ResourceLocator> genomeResources = GenomeManager.getInstance().getCurrentGenome().getAnnotationResources();
        Set<String> absoluteGenomeAnnotationPaths = genomeResources == null ? Collections.emptySet() :
                genomeResources.stream().map(rl -> rl.getPath()).collect(Collectors.toSet());

        if (dataFiles.size() > 0) {

            final List<String> errors = new ArrayList<String>();

            // Load files concurrently -- TODO, put a limit on # of threads?
            List<Thread> threads = new ArrayList(dataFiles.size());
            long t0 = System.currentTimeMillis();

            List<Runnable> synchronousLoads = new ArrayList<Runnable>();

            for (final ResourceLocator locator : dataFiles) {

                if (absoluteGenomeAnnotationPaths.contains(FileUtils.getAbsolutePath(locator.getPath(), genomePath))) {
                    continue;
                }

                Runnable runnable = () -> {
                    try {
                        // igv.load() loads and initializes tracks, but does not allocate them to panels.
                        List<Track> tracks = igv.load(locator);
                        for (Track track : tracks) {
                            if (track == null) {
                                log.info("Null track for resource " + locator.getPath());
                                continue;
                            }

                            String id = checkTrackId(track.getId());
                            if (id == null) {
                                log.info("Null track id for resource " + locator.getPath());
                                continue;
                            }

                            // id is often an absolute file path.  Use just the filename if unique
//                            if (!FileUtils.isRemote(id)) {
//                                String fn = (new File(id)).getName();
//                                if (!allTracks.containsKey(fn)) {
//                                    id = fn;
//                                }
//                            }

                            if (!allTracks.containsKey(id)) {
                                allTracks.put(id, new ArrayList<>());
                            }
                            allTracks.get(id).add(track);
                        }
                    } catch (Exception e) {
                        erroredResources.add(locator.getPath());
                        log.error("Error loading resource " + locator.getPath(), e);
                        String ms = "<b>" + locator.getPath() + "</b><br>&nbsp;&nbsp;" + e.toString() + "<br>";
                        errors.add(ms);
                    }
                };

                // Run synchronously if in batch mode or if there are no "track"  elements
                if (Globals.isBatch() || !hasTrackElments) {
                    synchronousLoads.add(runnable);
                } else {
                    Thread t = new Thread(runnable);
                    threads.add(t);
                    t.start();
                }
            }
            // Wait for all threads to complete
            for (Thread t : threads) {
                try {
                    t.join();
                } catch (InterruptedException ignore) {
                    log.error(ignore);
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
     */
    void processResource(Session session, Element element, String sessionPath) {

        String nodeName = element.getNodeName();
        boolean oldSession = nodeName.equals(SessionElement.DATA_FILE);

        String label = getAttribute(element, SessionAttribute.LABEL);
        String name = getAttribute(element, SessionAttribute.NAME);
        String sampleId = getAttribute(element, SessionAttribute.SAMPLE_ID);
        String description = getAttribute(element, SessionAttribute.DESCRIPTION);
        String type = getAttribute(element, SessionAttribute.TYPE);
        String coverage = getAttribute(element, SessionAttribute.COVERAGE);
        String mapping = getAttribute(element, SessionAttribute.MAPPING);
        String trackLine = getAttribute(element, SessionAttribute.TRACK_LINE);
        String colorString = getAttribute(element, SessionAttribute.COLOR);
        String index = getAttribute(element, SessionAttribute.INDEX);

        // URL to a database or webservice -- this is not common
        String serverURL = getAttribute(element, SessionAttribute.SERVER_URL);

        String path = getAttribute(element, SessionAttribute.PATH);

        // Older sessions used the "name" attribute for the path.
        if (oldSession && name != null) {
            path = name;
            int idx = name.lastIndexOf("/");
            if (idx > 0 && idx + 1 < name.length()) {
                name = name.substring(idx + 1);
            }
        }

        String absolutePath = (sessionPath == null || ParsingUtils.isDataURL(path)) ?
                path :
                getAbsolutePath(path, sessionPath);

        ResourceLocator resourceLocator = new ResourceLocator(serverURL, absolutePath);

        if (index != null) resourceLocator.setIndexPath(index);

        if (coverage != null) {
            String absoluteCoveragePath = coverage.equals(".") ? coverage : getAbsolutePath(coverage, sessionPath);
            resourceLocator.setCoverage(absoluteCoveragePath);
        }

        if (mapping != null) {
            String absoluteMappingPath = mapping.equals(".") ? mapping : getAbsolutePath(mapping, sessionPath);
            resourceLocator.setMappingPath(absoluteMappingPath);
        }

        String url = getAttribute(element, SessionAttribute.URL);
        if (url == null) {
            url = getAttribute(element, SessionAttribute.FEATURE_URL);
        }
        resourceLocator.setFeatureInfoURL(url);

        String infolink = getAttribute(element, SessionAttribute.HYPERLINK);
        if (infolink == null) {
            infolink = getAttribute(element, SessionAttribute.INFOLINK);
        }
        resourceLocator.setTrackInfoURL(infolink);


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
            resourceLocator.setFormat(type);
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

        process(session, elements, sessionPath);

    }

    private void processRegions(Session session, Element element, String sessionPath) {

        session.clearRegionsOfInterest();
        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);
    }

    private void processRegion(Session session, Element element, String sessionPath) {

        String chromosome = getAttribute(element, SessionAttribute.CHROMOSOME);
        String start = getAttribute(element, SessionAttribute.START_INDEX);
        String end = getAttribute(element, SessionAttribute.END_INDEX);
        String description = getAttribute(element, SessionAttribute.DESCRIPTION);

        RegionOfInterest region = new RegionOfInterest(chromosome, Integer.parseInt(start), Integer.parseInt(end), description);
        igv.addRegionOfInterest(region);

        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);
    }

    private void processHiddenAttributes(Session session, Element element) {
        NodeList elements = element.getChildNodes();
        Set<String> attributes = new HashSet();
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            if (childNode.getNodeName().equals("Attribute")) {
                attributes.add(((Element) childNode).getAttribute(SessionAttribute.NAME));
            }
        }
        session.setHiddenAttributes(attributes);
    }


    /**
     * For backward compatibility
     *
     * @param session
     * @param element
     */
    private void processVisibleAttributes(Session session, Element element) {

        NodeList elements = element.getChildNodes();
        if (elements.getLength() > 0) {
            Set<String> visibleAttributes = new HashSet();
            for (int i = 0; i < elements.getLength(); i++) {
                Node childNode = elements.item(i);
                if (childNode.getNodeName().equals(SessionElement.VISIBLE_ATTRIBUTE)) {
                    visibleAttributes.add(((Element) childNode).getAttribute(SessionAttribute.NAME));
                }
            }

            final List<String> attributeNames = AttributeManager.getInstance().getAttributeNames();
            Set<String> hiddenAttributes = new HashSet<String>(attributeNames);
            hiddenAttributes.removeAll(visibleAttributes);
            session.setHiddenAttributes(hiddenAttributes);

        }
    }

    private void processGeneList(Session session, Element element) {

        String name = getAttribute(element, SessionAttribute.NAME);

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
                if (childNode.getNodeName().equalsIgnoreCase(SessionElement.FRAME)) {
                    String frameName = getAttribute((Element) childNode, SessionAttribute.NAME);

                    ReferenceFrame f = frames.get(frameName);
                    if (f != null) {
                        reorderedFrames.add(f);
                        try {
                            String chr = getAttribute((Element) childNode, SessionAttribute.CHR);
                            final String startString =
                                    getAttribute((Element) childNode, SessionAttribute.START).replace(",", "");
                            final String endString =
                                    getAttribute((Element) childNode, SessionAttribute.END).replace(",", "");
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
        igv.resetFrames();
    }

    /**
     * Example <Filter match="all" name="" showTracks="false">
     *
     * @param session
     * @param element
     * @param sessionPath
     */
    private void processFilter(Session session, Element element, String sessionPath) {

        String match = getAttribute(element, SessionAttribute.FILTER_MATCH);
        String showAllTracks = getAttribute(element, SessionAttribute.FILTER_SHOW_ALL_TRACKS);

        String filterName = getAttribute(element, SessionAttribute.NAME);
        filter = new TrackFilter(filterName, null);

        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);

        // Save the filter
        session.setFilter(filter);

        // Update UI elements -- the matchAll and showAll state is kept in the UI, not the Filter object.  This seems wrong.  TODO
        if ("all".equalsIgnoreCase(match)) {
            igv.setFilterMatchAll(true);
        } else if ("any".equalsIgnoreCase(match)) {
            igv.setFilterMatchAll(false);
        }
        if ("true".equalsIgnoreCase(showAllTracks)) {
            igv.setFilterShowAllTracks(true);
        } else {
            igv.setFilterShowAllTracks(false);
        }
    }

    private void processFilterElement(Session session, Element element,
                                      String sessionPath) {

        if (filter == null) {
            throw new RuntimeException("Filter elements defined before filter");
        }

        String item = getAttribute(element, SessionAttribute.ITEM);
        String operator = getAttribute(element, SessionAttribute.OPERATOR);
        String value = getAttribute(element, SessionAttribute.VALUE);
        String booleanOperator = getAttribute(element, SessionAttribute.BOOLEAN_OPERATOR);

        Operator opEnum = CollUtils.findValueOf(Operator.class, operator);
        BooleanOperator boolEnum = BooleanOperator.valueOf(booleanOperator.toUpperCase());
        TrackFilterElement trackFilterElement = new TrackFilterElement(filter, item,
                opEnum, value, boolEnum);
        filter.add(trackFilterElement);

        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);
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
                (node.getNodeName().equalsIgnoreCase(SessionElement.DATA_TRACK) ||
                        node.getNodeName().equalsIgnoreCase(SessionElement.TRACK));
    }

    private void processPanel(Session session, Element element, String sessionPath) {

        if (panelElementPresent == false) {
            panelElementPresent = true;
            // First panel to be processed, do this only once.
            // Add any tracks loaded as a side effect of loading genome
            final List<Track> tmp = igv.getAllTracks();
            for (Track track : tmp) {
                String id = checkTrackId(track.getId());
                if (!allTracks.containsKey(id)) {
                    allTracks.put(id, new ArrayList<>());
                }
                allTracks.get(id).add(track);
            }

            // Tracks, if any,  will be placed back as prescribed by panel elements
            removeAllTracksFromPanels();
        }


        String panelName = element.getAttribute("name");
        if (panelName == null) {
            panelName = "Panel" + panelCounter++;
        }

        List<Track> panelTracks = new ArrayList();
        NodeList elements = element.getChildNodes();
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            if (nodeIsTrack(childNode)) {
                List<Track> tracks = processTrack(session, (Element) childNode, sessionPath);
                if (tracks != null) {
                    panelTracks.addAll(tracks);
                }
            } else {
                process(session, childNode, sessionPath);
            }
        }

        //We make a second pass through, resolving references to tracks which may have been processed after being referenced.
        //For instance if Track 2 referenced Track 4
        for (Track track : panelTracks) {
            if (track instanceof FeatureTrack) {
                FeatureTrack featureTrack = (FeatureTrack) track;
                featureTrack.updateTrackReferences(panelTracks);
            } else if (track instanceof DataSourceTrack) {
                DataSourceTrack dataTrack = (DataSourceTrack) track;
                dataTrack.updateTrackReferences(panelTracks);
            }
        }

        TrackPanel panel = igv.getTrackPanel(panelName);
        panel.addTracks(panelTracks);
    }

    private void processPanelLayout(Session session, Element element) {

        String nodeName = element.getNodeName();

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
     * @return
     */

    private List<Track> processTrack(Session session, Element element, String sessionPath) {

        String id = checkTrackId(getAttribute(element, SessionAttribute.ID));

        // Find track matching element id, created earlier from "Resource or File" elements, or during genome load.
        // Normally this is a single track, but that can't be assumed as uniqueness of "id" is not enforced.
        List<Track> matchedTracks = getTracksById(id);

        if (matchedTracks == null) {
            //Try creating an absolute path for the id
            if (id != null) {
                matchedTracks = allTracks.get(getAbsolutePath(id, sessionPath));
            }
        }

        if (matchedTracks != null) {
            for (final Track track : matchedTracks) {
                track.unmarshalXML(element, version);
                allocatedToPanel.add(track);
            }
        } else {
            // No match found, element represents either (1) a track from a file that errored during load, for example
            // a file that has been deleted, or (2) a track not created from "Resource" or genome load.  These include
            // reference sequence, combined,  and merged tracks.

            String className = getAttribute(element, "clazz");
            if (resourceIndpendentTracks.stream().anyMatch(c -> className.contains(c))) {
                try {
                    Track track = createTrack(className, element);
                    if (track != null) {

                        track.unmarshalXML(element, version);
                        matchedTracks = Arrays.asList(track);
                        allTracks.put(checkTrackId(track.getId()), matchedTracks);   // Important for second pass

                        // Special tracks
                        if (className.contains("CombinedDataTrack")) {
                            combinedDataSourceTracks.add(new Pair(track, element));
                        }

                        if (className.contains("MergedTracks")) {
                            List<DataTrack> memberTracks = processChildTracks(session, element, sessionPath);
                            ((MergedTracks) track).setMemberTracks(memberTracks);

                            NodeList nodeList = element.getElementsByTagName("DataRange");
                            if (nodeList != null && nodeList.getLength() > 0) {
                                Element dataRangeElement = (Element) nodeList.item(0);
                                try {
                                    DataRange dataRange = new DataRange(dataRangeElement, version);
                                    track.setDataRange(dataRange);
                                } catch (Exception e) {
                                    log.error("Unrecognized DataRange");
                                }
                            }
                        }
                    } else {
                        log.warn("Warning.  No tracks were found with id: " + id + " in session file");
                    }
                } catch (Exception e) {
                    log.error("Error restoring track: " + element.toString(), e);
                    MessageUtils.showMessage("Error loading track: " + element.toString());
                }
            }
        }

        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);

        return matchedTracks;
    }

    /**
     * Recursively loop through children of {@code element}, and process them as tracks iff they are determined
     * to be so
     *
     * @param session
     * @param element
     * @param sessionPath
     * @return List of processed tracks.
     */
    private List<DataTrack> processChildTracks(Session session, Element element, String sessionPath) {

        NodeList memberTrackNodes = element.getChildNodes();
        List<DataTrack> memberTracks = new ArrayList<>(memberTrackNodes.getLength());
        for (int index = 0; index < memberTrackNodes.getLength(); index++) {
            Node memberNode = memberTrackNodes.item(index);
            if (nodeIsTrack(memberNode)) {
                List<Track> addedTracks = processTrack(session, (Element) memberNode, sessionPath);
                if (addedTracks != null) {
                    for (Track t : addedTracks) {
                        if (t instanceof DataTrack) {
                            memberTracks.add((DataTrack) t);
                        } else {
                            log.error("Unexpected MergedTrack member class: " + t.getClass().getName());
                        }
                    }
                }
            }
        }
        return memberTracks;
    }

    /**
     * Process combined data tracks, these are tracks composed from dependent tracks by simple arithmetic operations.
     */
    private void processCombinedDataSourceTracks() {

        Map<CombinedDataTrack, CombinedDataSource> sourceMap = new HashMap<>();

        // First pass -- create data sources
        for (Pair<CombinedDataTrack, Element> pair : this.combinedDataSourceTracks) {

            Element element = pair.getSecond();
            CombinedDataTrack combinedTrack = pair.getFirst();

            DataTrack track1 = null;
            DataTrack track2 = null;
            List<Track> tmp = getTracksById(element.getAttribute("track1"));
            if (tmp != null && tmp.size() > 0) {
                track1 = (DataTrack) tmp.get(0);
            }
            tmp = getTracksById(element.getAttribute("track2"));
            if (tmp != null && tmp.size() > 0) {
                track2 = (DataTrack) tmp.get(0);
            }

            if (track1 == null || track2 == null) {
                log.error("Missing track for combined track: " + pair.getFirst().getName());
                return;
            }
            CombinedDataSource.Operation op = CombinedDataSource.Operation.valueOf(element.getAttribute("op"));

            CombinedDataSource source = new CombinedDataSource(track1, track2, op);
            sourceMap.put(combinedTrack, source);
        }

        // Now set datasource on tracks.  This needs to be deferred as combined data sources can reference other
        // combined sources, we need to instantiate the entire tree before using a datasource
        for (Map.Entry<CombinedDataTrack, CombinedDataSource> entry : sourceMap.entrySet()) {
            entry.getKey().setDatasource(entry.getValue());
        }

    }

    private void processColorScales(Session session, Element element, String sessionPath) {

        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);
    }

    private void processColorScale(Session session, Element element, String sessionPath) {

        String trackType = getAttribute(element, SessionAttribute.TYPE);
        String value = getAttribute(element, SessionAttribute.VALUE);

        setColorScaleSet(session, trackType, value);

        NodeList elements = element.getChildNodes();
        process(session, elements, sessionPath);
    }

    private void processPreferences(Session session, Element element) {

        NodeList elements = element.getChildNodes();
        for (int i = 0; i < elements.getLength(); i++) {
            Node child = elements.item(i);
            if (child.getNodeName().equalsIgnoreCase(SessionElement.PROPERTY)) {
                Element childNode = (Element) child;
                String name = getAttribute(childNode, SessionAttribute.NAME);
                String value = getAttribute(childNode, SessionAttribute.VALUE);
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
    private void process(Session session, NodeList elements, String sessionPath) {
        for (int i = 0; i < elements.getLength(); i++) {
            Node childNode = elements.item(i);
            process(session, childNode, sessionPath);
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


    /**
     * Uses {@link #currentReader} to lookup matching tracks by id, or
     * searches allTracks if sessionReader is null
     *
     * @param trackId
     * @param allTracks
     * @return
     */
    public static Track getMatchingTrack(String trackId, List<Track> allTracks) {
        trackId = checkTrackId(trackId);
        IGVSessionReader reader = currentReader.get();
        List<Track> matchingTracks;
        if (reader != null) {
            matchingTracks = reader.getTracksById(trackId);
        } else {
            if (allTracks == null)
                throw new IllegalStateException("No session reader and no tracks to search to resolve Track references");
            matchingTracks = new ArrayList<>();
            for (Track track : allTracks) {
                if (trackId.equals(checkTrackId(track.getId()))) {
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

    /**
     * Create a track object for the given class name.  In the past this was done by reflection.   This was a bad
     * idea.  Among other issues the full class name gets hardcoded into sessions preventing future refactoring
     * or renaming.  Thus the lookup table approach below, with reflection at the end in case we miss any.
     *
     * @param className
     * @return
     * @throws ClassNotFoundException
     * @throws InstantiationException
     * @throws IllegalAccessException
     * @throws java.lang.reflect.InvocationTargetException
     * @throws NoSuchMethodException
     */
    private Track createTrack(String className, Element element) throws ClassNotFoundException, InstantiationException, IllegalAccessException, java.lang.reflect.InvocationTargetException, NoSuchMethodException {

        if (className.contains("BasePairTrack")) {
            return new BasePairTrack();
        } else if (className.contains("BlatTrack")) {
            return new BlatTrack();
        } else if (className.contains("ClusterTrack")) {
            return new ClusterTrack();
        } else if (className.contains("CNFreqTrack")) {
            return new CNFreqTrack();
        } else if (className.contains("CoverageTrack")) {
            return new CoverageTrack();
        } else if (className.contains("DataSourceTrack")) {
            return new DataSourceTrack();
        } else if (className.contains("DSITrack")) {
            return new DSITrack();
        } else if (className.contains("EWigTrack")) {
            return new EWigTrack();
        } else if (className.contains("FeatureTrack")) {
            return new FeatureTrack();
        } else if (className.contains("MotifTrack")) {
            return new MotifTrack();
        } else if (className.contains("GisticTrack")) {
            return new GisticTrack();
        } else if (className.contains("InteractionTrack")) {
            return new InteractionTrack();
        } else if (className.contains("MergedTracks")) {
            return new MergedTracks();
        } else if (className.contains("CombinedDataTrack")) {
            String id = element.getAttribute("id");
            String name = element.getAttribute("name");
            return new CombinedDataTrack(id, name);
        } else if (className.contains("MultipleAlignmentTrack")) {
            return new MultipleAlignmentTrack();
        } else if (className.contains("MutationTrack")) {
            return new MutationTrack();
        } else if (className.contains("SelectableFeatureTrack")) {
            return new SelectableFeatureTrack();
        } else if (className.contains("SpliceJunctionTrack")) {
            return new SpliceJunctionTrack();
        } else if (className.contains("VariantTrack")) {
            return new VariantTrack();
        } else if (className.contains("SequenceTrack")) {
            return new SequenceTrack("Reference sequence");
        } else {
            log.warn("Unrecognized class name: " + className);
            try {
                Class clazz = SessionElement.getClass(className);
                return (Track) clazz.getConstructor().newInstance();
            } catch (Exception e) {
                log.error("Error attempting Track creation ", e);
                return null;
            }
        }
    }

    private void removeAllTracksFromPanels() {
        List<TrackPanel> panels = igv.getTrackPanels();
        for (TrackPanel trackPanel : panels) {
            trackPanel.removeAllTracks();
        }
    }

    private String getAbsolutePath(String path, String sessionPath) {
        String absolute = FileUtils.getAbsolutePath(path, sessionPath);
        if (!(new File(absolute)).exists() && path.startsWith("~")) {
            String tmp = FileUtils.getAbsolutePath(path.replaceFirst("~", System.getProperty("user.home")), sessionPath);
            if ((new File(tmp)).exists()) {
                absolute = tmp;
            }
        }
        return absolute;
    }

    public List<Track> getTracksById(String trackId) {
        List<Track> tracks = allTracks.get(trackId);
        if (tracks == null) {
            // ID is usually a full file path or URL.  See if we can find the track with just filename
            if (!FileUtils.isRemote(trackId)) {
                String fn = (new File(trackId)).getName();
                tracks = allTracks.get(fn);
            }
            if (tracks == null) {
                // If still no match search for legacy gene annotation specifier
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                String legacyGeneTrackID = genome.getId() + "_genes";
                if (trackId.equals(legacyGeneTrackID)) {
                    if (genome.getGeneTrack() != null) {
                        return Arrays.asList(genome.getGeneTrack());
                    }
                }
            }
        }
        return tracks;
    }

    private static String checkTrackId(String id) {
        return id == null ? null : (geneTrackIds.containsKey(id)) ? geneTrackIds.get(id) : id;
    }

    /**
     * Set of track classes that are not backed by resources (files).
     */
    private static Set<String> resourceIndpendentTracks = new HashSet<>(Arrays.asList("MergedTracks", "CombinedDataTrack",
            "SequenceTrack", "BlatTrack", "MotifTrack"));

    /**
     * Some synonyms for gene tracks with multiple URLs over time.  Allows matching old sessions to updated genomes.
     */
    private static Map<String, String> geneTrackIds = Map.of(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz", "hg38_gene_track",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz", "hg38_gene_track",
            "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.txt.gz", "hg38_gene_track",
            "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz", "hg38_gene_track",
            "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefGene.txt.gz", "hg38_gene_track",
            "https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg38/ncbiRefGene.txt.gz", "hg38_gene_track",
            "https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz", "hg19_gene_track",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz", "hg19_gene_track");
}
