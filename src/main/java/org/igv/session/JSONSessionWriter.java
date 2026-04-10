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

        // Gene list (multi-locus view)
        if (FrameManager.isGeneListMode()) {
            GeneList geneList = session.getCurrentGeneList();
            if (geneList != null) {
                JSONObject geneListJson = new JSONObject();
                geneListJson.put("name", geneList.getName());

                JSONArray lociArray = new JSONArray();
                for (String l : geneList.getLoci()) {
                    lociArray.put(l);
                }
                geneListJson.put("loci", lociArray);

                // Store the current frame extents, which may differ from the original gene list loci
                // (e.g. due to zooming or panning)
                JSONArray framesArray = new JSONArray();
                for (ReferenceFrame frame : FrameManager.getFrames()) {
                    JSONObject frameJson = new JSONObject();
                    frameJson.put("name", frame.getName());
                    frameJson.put("chr", frame.getChrName());
                    frameJson.put("start", frame.getOrigin());
                    frameJson.put("end", frame.getEnd());
                    framesArray.put(frameJson);
                }
                geneListJson.put("frames", framesArray);

                sessionObject.put("geneList", geneListJson);
            }
        }

        // Hidden attributes
        Set<String> hiddenAttributes = session.getHiddenAttributes();
        if (hiddenAttributes != null && !hiddenAttributes.isEmpty()) {
            JSONArray hiddenArray = new JSONArray();
            for (String attr : hiddenAttributes) {
                hiddenArray.put(attr);
            }
            sessionObject.put("hiddenAttributes", hiddenArray);
        }

        return sessionObject.toString(2);

    }




    private boolean isUseRelative(File outputFile) {
        return outputFile != null &&
                PreferencesManager.getPreferences().getAsBoolean(Constants.SESSION_RELATIVE_PATH);
    }

}

