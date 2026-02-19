package org.igv.session;

import org.igv.feature.RegionOfInterest;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.lists.GeneList;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.panel.TrackPanel;
import org.igv.util.*;
import org.json.JSONArray;
import org.json.JSONObject;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author jrobinso
 */
public class JSONSessionWriter {

    static Logger log = LogManager.getLogger(JSONSessionWriter.class);

    private IGV igv;
    private Session session;
    private static int CURRENT_VERSION = 8;
    private File outputFile;
    private Document document;

    public JSONSessionWriter(IGV igv) {
        this.igv = igv;
    }

    /**
     * Save the session as a JSON file
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
            log.error("Session Management Error: NULL outputFile");
        }

        String json = createJsonFromSession(session);

        Writer fileWriter = null;
        try {
            fileWriter = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(outputFile), "UTF8"));
            fileWriter.write(json);
        } finally {
            if (fileWriter != null) {
                fileWriter.close();
            }
        }
    }

    public String createJsonFromSession(Session session) throws RuntimeException {

        this.session = session;

        JSONObject sessionObject = new JSONObject();

        JSONObject genomeJson = GenomeManager.getInstance().getCurrentGenome().getConfig().toJSON();
        genomeJson.remove("tracks");   // Don't include tracks in the genome definition, they are included in the session
        sessionObject.put("reference", genomeJson);

        String locus = session.getCurrentLocus();
        if (locus != null && !FrameManager.isGeneListMode()) {
            sessionObject.put(SessionAttribute.LOCUS, locus);
        }

        String groupBy = session.getGroupByAttribute();
        if (groupBy != null) {
            sessionObject.put(SessionAttribute.GROUP_TRACKS_BY, groupBy);
        }

        int nextAutoscaleGroup = session.getNextAutoscaleGroup();
        if (nextAutoscaleGroup > 1) {
            sessionObject.put(SessionAttribute.NEXT_AUTOSCALE_GROUP, String.valueOf(nextAutoscaleGroup));
        }

        JSONArray tracks = new JSONArray();
        for (Track track : IGV.getInstance().getAllTracks()) {

            JSONObject trackJson = new JSONObject();

            track.marshalJSON(trackJson);

            tracks.put(trackJson);
        }

        sessionObject.put("tracks", tracks);

        return sessionObject.toString(2);

    }

    private void writeFilters(Session session, Element globalElement, Document document) {
        TrackFilter trackFilter = session.getFilter();
        if (trackFilter != null) {

            Element filter = document.createElement(SessionElement.FILTER);

            if (!IGV.getInstance().isFilterMatchAll()) {
                filter.setAttribute(SessionAttribute.FILTER_MATCH, "any");
            } else {    // Defaults to match all
                filter.setAttribute(SessionAttribute.FILTER_MATCH, "all");
            }

            if (IGV.getInstance().isFilterShowAllTracks()) {
                filter.setAttribute(SessionAttribute.FILTER_SHOW_ALL_TRACKS, "true");
            } else {    // Defaults
                filter.setAttribute(SessionAttribute.FILTER_SHOW_ALL_TRACKS, "false");
            }
            globalElement.appendChild(filter);

            // Process FilterElement elements
            Iterator iterator = session.getFilter().getFilterElements();
            while (iterator.hasNext()) {

                FilterElement trackFilterElement = (FilterElement) iterator.next();

                Element filterElementElement =
                        document.createElement(SessionElement.FILTER_ELEMENT);
                filterElementElement.setAttribute(SessionAttribute.ITEM,
                        trackFilterElement.getAttributeKey());
                filterElementElement.setAttribute(
                        SessionAttribute.OPERATOR,
                        trackFilterElement.getComparisonOperator().getValue());
                filterElementElement.setAttribute(SessionAttribute.VALUE,
                        trackFilterElement.getValue());

                filter.appendChild(filterElementElement);
            }
        }
    }

    private void writeRegionsOfInterest(Element globalElement, Document document) {
        Collection<RegionOfInterest> regions = session.getAllRegionsOfInterest();
        if ((regions != null) && !regions.isEmpty()) {

            Element regionsElement = document.createElement(SessionElement.REGIONS);
            for (RegionOfInterest region : regions) {
                Element regionElement = document.createElement(SessionElement.REGION);
                regionElement.setAttribute(SessionAttribute.CHROMOSOME, region.getChr());
                regionElement.setAttribute(SessionAttribute.START_INDEX, String.valueOf(region.getStart()));
                regionElement.setAttribute(SessionAttribute.END_INDEX, String.valueOf(region.getEnd()));
                if (region.getDescription() != null) {
                    regionElement.setAttribute(SessionAttribute.DESCRIPTION, region.getDescription());
                }
                regionsElement.appendChild(regionElement);
            }
            globalElement.appendChild(regionsElement);
        }
    }

    private void writeHiddenAttributes(Session session, Element globalElement, Document document) {
        Element hiddenAttributes = document.createElement(SessionElement.HIDDEN_ATTRIBUTES);
        for (String attribute : session.getHiddenAttributes()) {
            Element regionElement = document.createElement(SessionElement.ATTRIBUTE);
            regionElement.setAttribute(SessionAttribute.NAME, attribute);
            hiddenAttributes.appendChild(regionElement);
        }
        globalElement.appendChild(hiddenAttributes);
    }

    private void writeGeneList(Element globalElement, Document document) {

        GeneList geneList = session.getCurrentGeneList();

        if (geneList != null) {

            Element geneListElement = document.createElement(SessionElement.GENE_LIST);
            geneListElement.setAttribute(SessionAttribute.NAME, geneList.getName());
            StringBuffer genes = new StringBuffer();
            for (String gene : geneList.getLoci()) {
                genes.append(gene);
                genes.append("\n");
            }
            geneListElement.setTextContent(genes.toString());
            globalElement.appendChild(geneListElement);


            // Now store the list of frames visible.  This seems redundant, but frame extent can be changed after
            // "gene list" definition, for example by zooming out or panning in a frame
            for (ReferenceFrame frame : FrameManager.getFrames()) {
                Element frameElement = document.createElement(SessionElement.FRAME);
                frameElement.setAttribute(SessionAttribute.NAME, frame.getName());
                frameElement.setAttribute(SessionAttribute.CHR, frame.getChrName());
                frameElement.setAttribute(SessionAttribute.START, String.valueOf(frame.getOrigin()));
                frameElement.setAttribute(SessionAttribute.END, String.valueOf(frame.getEnd()));

                geneListElement.appendChild(frameElement);

            }
        }
    }

    private boolean isUseRelative(File outputFile) {
        return outputFile != null &&
                PreferencesManager.getPreferences().getAsBoolean(Constants.SESSION_RELATIVE_PATH);
    }

    private void writePanels(Element globalElement, Document document) throws DOMException {

        for (TrackPanel trackPanel : IGV.getInstance().getTrackPanels()) {

            // TODO -- loop through panels groups, rather than skipping groups to tracks

            List<Track> tracks = trackPanel.getTracks();
            if ((tracks != null) && !tracks.isEmpty()) {

                Element panelElement = document.createElement(SessionElement.PANEL);
                panelElement.setAttribute("name", trackPanel.getName());
                panelElement.setAttribute("height", String.valueOf(trackPanel.getHeight()));
                panelElement.setAttribute("width", String.valueOf(trackPanel.getWidth()));

                for (Track track : tracks) {

                    Element element = document.createElement("Track");
                    element.setAttribute("clazz", SessionElement.getXMLClassName(track.getClass()));

                    String id = track.getId();
                    if (isUseRelative(outputFile) && !FileUtils.isRemote(id)) {
                        id = FileUtils.getRelativePath(outputFile.getAbsolutePath(), id);
                    }
                    element.setAttribute("id", id);

                    track.marshalXML(document, element);

                    if (track.isNumeric() && track.getDataRange() != null) {
                        Element dataRangeElement = document.createElement(SessionElement.DATA_RANGE);
                        track.getDataRange().marshalXML(document, dataRangeElement);
                        element.appendChild(dataRangeElement);
                    }

                    panelElement.appendChild(element);

                }

                globalElement.appendChild(panelElement);
            }
        }
    }

    private void writePanelLayout(Element globalElement, Document document) {

    }

    /**
     * @return A set of the load data files.
     */
    public Collection<ResourceLocator> getResourceLocatorSet() {

        Collection<ResourceLocator> locators = new ArrayList();

        Collection<ResourceLocator> currentTrackFileLocators =
                IGV.getInstance().getDataResourceLocators();

        if (currentTrackFileLocators != null) {

            // Filter data files that are included in genome annotations
            final Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
            if (currentGenome != null) {
                List<ResourceLocator> genomeResources = currentGenome.getAnnotationResources();
                Set<String> absoluteGenomeAnnotationPaths = genomeResources == null ? Collections.emptySet() :
                        genomeResources.stream().map(rl -> rl.getPath()).collect(Collectors.toSet());

                for (ResourceLocator locator : currentTrackFileLocators) {
                    if (!absoluteGenomeAnnotationPaths.contains(locator.getPath())) {
                        locators.add(locator);
                    }
                }
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

