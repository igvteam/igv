package org.igv.ui.action;

import org.igv.encode.EncodeTrackChooserFactory;
import org.igv.encode.FileRecord;
import org.igv.encode.TrackChooser;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.ui.WaitCursorManager;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * Date: 11/2/13
 * Time: 6:39 PM
 */
public class BrowseEncodeAction extends MenuAction {


    public enum Type {
        UCSC,
        SIGNALS_CHIP,
        SIGNALS_OTHER,
        OTHER,
        HIC,
        FOUR_DN
    }

    private static Logger log = LogManager.getLogger(BrowseEncodeAction.class);

    private static Map<String, Color> colors;

    static {
        colors = new HashMap<>();
        colors.put("H3K27AC", new Color(200, 0, 0));
        colors.put("H3K27ME3", new Color(200, 0, 0));
        colors.put("H3K36ME3", new Color(0, 0, 150));
        colors.put("H3K4ME1", new Color(0, 150, 0));
        colors.put("H3K4ME2", new Color(0, 150, 0));
        colors.put("H3K4ME3", new Color(0, 150, 0));
        colors.put("H3K9AC", new Color(100, 0, 0));
        colors.put("H3K9ME1", new Color(100, 0, 0));
    }

    /** Sample info attributes
     *      * Properties available in various sets
     *      *   UCSC Encode:  path	cell	dataType	antibody	view	replicate	type	lab	hub
     *      *   Encode:       Biosample	AssayType	Target	BioRep	TechRep	OutputType
     *      *   4DN:          Type	Biosource	Assay	Replicate	Dataset Accession	Experiment	name
     */
    static Set<String> sampleInfoAttributes = new HashSet<>(Arrays.asList(
            "dataType", "cell", "antibody", "lab", "Biosample", "AssayType", "Target", "Biosource"));

    static Set<String> trackLineAttributes = new HashSet<>(Arrays.asList(
            "name", "description", "color", "altColor", "visibility", "maxHeightPixels",
            "viewLimits", "autoScale", "priority"));

    private final Type type;

    IGV igv;

    public BrowseEncodeAction(String label, int mnemonic, Type type, IGV igv) {
        super(label, null, mnemonic);
        this.type = type;
        this.igv = igv;
    }


    @Override
    public void actionPerformed(ActionEvent event) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        final WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
        SwingWorker worker = new SwingWorker<TrackChooser, Void>() {
            @Override
            protected TrackChooser doInBackground() throws Exception {
                return EncodeTrackChooserFactory.getInstance(genome.getId(), BrowseEncodeAction.this.type);
            }

            @Override
            protected void done() {
                WaitCursorManager.removeWaitCursor(token);
                try {
                    TrackChooser chooser = get();
                    if (chooser == null) {
                        MessageUtils.showMessage("Encode data is not available for " + genome.getDisplayName() + " through IGV.");
                        return;
                    }

                    chooser.setVisible(true);
                    if (chooser.isCanceled()) return;


                    java.util.List<FileRecord> records = chooser.getAllRecords();

                    final List<Track> loadedTracks = IGV.getInstance().getAllTracks().stream().filter(t -> t.getResourceLocator() != null).toList();
                    final Set<String> loadedTrackPaths = new HashSet<>(loadedTracks.stream().map(t -> t.getResourceLocator().getPath()).toList());
                    final Set<String> trackPathsToRemove = new HashSet<>();
                    final List<ResourceLocator> tracksToLoad = new ArrayList<>(records.size());

                    for (FileRecord record : records) {
                        if (record.isSelected()) {
                            if (!loadedTrackPaths.contains(record.getPath())) {
                                tracksToLoad.add(getResourceLocator(record));
                            }
                        } else {
                            trackPathsToRemove.add(record.getPath());
                        }
                    }

                    List<Track> tracksToRemove = loadedTracks.stream().filter(t -> trackPathsToRemove.contains(t.getResourceLocator().getPath())).toList();
                    igv.deleteTracks(tracksToRemove);

                    igv.loadTracks(tracksToLoad);

                } catch (Exception e) {
                    log.error("Error opening Encode browser", e);
                    throw new RuntimeException(e);
                }
            }
        };
        worker.execute();
    }

    private ResourceLocator getResourceLocator(FileRecord record) {
        ResourceLocator rl = new ResourceLocator(record.getPath());
        rl.setName(getTrackName(record));
        Map<String, String> attributes = record.getAttributes();
        StringBuffer trackLine = new StringBuffer();

        String antibody = attributes.containsKey("antibody") ? attributes.get("antibody") : attributes.get("Target");
        if (antibody != null) {
            rl.setColor(colors.get(antibody.toUpperCase()));
        }

        for (Map.Entry<String, String> entry : attributes.entrySet()) {
            String value = entry.getValue();
            String normalizedKey = normalizeAttributeName(entry.getKey());
            if (value != null && value.length() > 0 && sampleInfoAttributes.contains(normalizedKey)) {
                AttributeManager.getInstance().addAttribute(rl.getName(),normalizedKey, value);
            } else if (value != null && value.length() > 0 && trackLineAttributes.contains(normalizedKey)) {
                trackLine.append(normalizedKey + "=\"" + value + "\" ");
            }
        }
        rl.setMetadata(attributes);
        if (trackLine.length() > 0) {
            rl.setTrackLine(trackLine.toString().trim());
        }

        return rl;
    }

    /**
     * Normalize attribute names across different data sources (Encode, UCSC, 4DN)
     * @param rawName
     * @return
     */
    private static String normalizeAttributeName(String rawName) {
        switch (rawName) {
            case "cell":
            case "Biosource":
                return "Biosample";
            case "antibody":
                return "Target";
            default:
                return rawName;
        }
    }


    /**
     * Return a friendly name for the track.
     *
     * Properties available in various sets
     *   UCSC Encode:  path	cell	dataType	antibody	view	replicate	type	lab	hub
     *   Encode:       Biosample	AssayType	Target	BioRep	TechRep	OutputType
     *   4DN:          Type	Biosource	Assay	Replicate	Dataset Accession	Experiment	name
     *
     * @return
     */
    public String getTrackName(FileRecord record) {

        Map<String, String> attributes = record.getAttributes();
        if (attributes.containsKey("name")) {
            return attributes.get("name");
        }
        if(attributes.containsKey("Dataset")) {
            return attributes.get("Dataset");
        }

        StringBuffer sb = new StringBuffer();
        if (attributes.containsKey("cell")) sb.append(attributes.get("cell"));
        else if (attributes.containsKey("Biosample")) sb.append(attributes.get("Biosample"));

        if (attributes.containsKey("antibody")) sb.append(" " + attributes.get("antibody"));
        else if (attributes.containsKey("Target")) sb.append(" " + attributes.get("Target"));

        if (attributes.containsKey("dataType")) sb.append(" " + attributes.get("dataType"));
        else if (attributes.containsKey("AssayType")) sb.append(" " + attributes.get("AssayType"));
        else if (attributes.containsKey("Assay")) sb.append(" " + attributes.get("Assay"));

        if (attributes.containsKey("view")) sb.append(" " + attributes.get("view"));
        else if (attributes.containsKey("OutputType")) sb.append(" " + attributes.get("OutputType"));

        String trackName = sb.toString().trim();
        if (sb.length() == 0) {
            trackName = (new File(record.getPath())).getName();
        }

        return trackName;

    }
}
