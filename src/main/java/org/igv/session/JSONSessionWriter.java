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
            sessionObject.put(SessionAttribute.GROUP_SAMPLES_ATTRIBUTE, groupBy);
        }

        int nextAutoscaleGroup = session.getNextAutoscaleGroup();
        if (nextAutoscaleGroup > 1) {
            sessionObject.put(SessionAttribute.NEXT_AUTOSCALE_GROUP, String.valueOf(nextAutoscaleGroup));
        }

        Set<ResourceLocator> sampleInfoLocators = AttributeManager.getInstance().getLoadedResources();
        if(sampleInfoLocators != null && !sampleInfoLocators.isEmpty()) {
            JSONArray sampleInfoArray = new JSONArray();
            for(ResourceLocator locator : sampleInfoLocators) {
                JSONObject sampleInfoJson = new JSONObject();
                sampleInfoJson.put("name", locator.getName());
                sampleInfoJson.put("url", locator.getPath());
                sampleInfoJson.put("format", "sampleInfo");
                sampleInfoArray.put(sampleInfoJson);
            }
            sessionObject.put("sampleinfo", sampleInfoArray);
        }

        Collection<RegionOfInterest> regions = session.getAllRegionsOfInterest();
        if ((regions != null) && !regions.isEmpty()) {
            JSONArray roiArray = new JSONArray();
            JSONObject roiObject = new JSONObject();
            JSONArray featuresArray = new JSONArray();

            for (RegionOfInterest region : regions) {
                JSONObject featureJson = new JSONObject();
                featureJson.put("chr", region.getChr());
                featureJson.put("start", region.getStart());
                featureJson.put("end", region.getEnd());
                if (region.getDescription() != null) {
                    featureJson.put("description", region.getDescription());
                }
                featuresArray.put(featureJson);
            }

            roiObject.put("features", featuresArray);
            roiArray.put(roiObject);
            sessionObject.put("roi", roiArray);
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

