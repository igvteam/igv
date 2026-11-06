package org.igv.figwig;

import org.igv.Globals;
import org.igv.encode.FileRecord;
import org.igv.encode.TrackChooser;
import org.igv.ui.IGV;
import org.igv.ui.util.UIUtilities;

import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.function.Function;

/**
 * Which tracks of a figwig file to load.  A file of many tracks (the cells of a single-cell
 * experiment, every experiment of a cell type) opens IGV's track chooser, the table the ENCODE
 * browser uses, with a row per track and a column per attribute.  Smaller files, headless and batch
 * runs, and session restores, which name their own tracks, load every track.
 */
public class FigWigTrackChooser {

    /** Files with more tracks than this open the chooser. */
    public static final int CHOOSER_MIN_TRACKS = 50;

    /** For tests: replaces the dialog, and is used whatever the file's size or the mode. */
    public static Function<FigWigFile, Set<String>> override = null;

    /**
     * The names of the tracks to load, or null for every track.  An empty set means the user
     * cancelled.
     */
    public static Set<String> choose(FigWigFile file, String title) {
        if (override != null) {
            return override.apply(file);
        }
        if (file.nTracks() <= CHOOSER_MIN_TRACKS || Globals.isHeadless() || Globals.isBatch() || !IGV.hasInstance()
                || IGV.getInstance().isRestoringSession()) {
            return null;
        }
        List<String> headings = headings(file);
        List<FileRecord> records = records(file);
        Set<String> chosen = new LinkedHashSet<>();
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            Frame owner = IGV.getInstance().getMainFrame();
            TrackChooser chooser = new TrackChooser(owner, headings, records, title);
            Rectangle bounds = owner.getBounds();
            chooser.setSize(new Dimension(Math.min(bounds.width, 1200), Math.min(bounds.height, 800)));
            chooser.setLocationRelativeTo(owner);
            chooser.setVisible(true);
            if (!chooser.isCanceled()) {
                for (FileRecord r : chooser.getSelectedRecords()) {
                    chosen.add(r.getAttributeValue("name"));
                }
            }
        });
        return chosen;
    }

    /** The chooser's columns: the track's name, then its attributes in order of first appearance, at most 15. */
    static List<String> headings(FigWigFile file) {
        LinkedHashSet<String> names = new LinkedHashSet<>();
        names.add("name");
        for (Map<String, String> attributes : file.trackAttributes) {
            names.addAll(attributes.keySet());
        }
        List<String> headings = new ArrayList<>(names);
        return headings.size() > 15 ? headings.subList(0, 15) : headings;
    }

    static List<FileRecord> records(FigWigFile file) {
        List<FileRecord> records = new ArrayList<>(file.nTracks());
        for (int t = 0; t < file.nTracks(); t++) {
            Map<String, String> attributes = new LinkedHashMap<>();
            attributes.put("name", file.trackNames[t]);
            attributes.putAll(file.trackAttributes.get(t));
            records.add(new FileRecord(file.trackNames[t], attributes));
        }
        return records;
    }
}
