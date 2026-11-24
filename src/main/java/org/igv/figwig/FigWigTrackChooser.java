package org.igv.figwig;

import org.figwig.FigWigFile;
import org.igv.Globals;
import org.igv.encode.FileRecord;
import org.igv.encode.TrackChooser;
import org.igv.ui.IGV;
import org.igv.ui.util.UIUtilities;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * Which tracks of a figwig file to load.  A file of many tracks (the cells of a single-cell
 * experiment, every experiment of a cell type) opens IGV's track chooser, the table the ENCODE
 * browser uses, with a row per track and a column per attribute.  It is asked for by the user
 * interface before a file is loaded (IGV.loadTracks), and the choice goes to the loader on the
 * ResourceLocator.  Smaller files and batch runs load every track; so do session restores, which do
 * not come this way.
 */
public class FigWigTrackChooser {

    /** Files with more tracks than this open the chooser. */
    public static final int CHOOSER_MIN_TRACKS = 50;

    /**
     * Ask which tracks of a figwig file with many tracks to load, and set them on the locator.  Does
     * nothing for other files.  False when the user cancels, and the file is then not to be loaded.
     */
    public static boolean choose(ResourceLocator locator) throws IOException {
        String format = locator.getFormat();
        if (!("fw".equals(format) || "figwig".equals(format)) || Globals.isBatch()) {
            return true;
        }
        List<String> headings;
        List<FileRecord> records;
        FigWigFile file = IGVFigWigByteSource.open(locator.getPath());
        try {
            if (file.nTracks() <= CHOOSER_MIN_TRACKS) {
                return true;
            }
            headings = headings(file);
            records = records(file);
        } finally {
            file.close();
        }
        Set<String> chosen = new LinkedHashSet<>();
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            Frame owner = IGV.getInstance().getMainFrame();
            TrackChooser chooser = new TrackChooser(owner, headings, records, locator.getTrackName());
            JLabel header = new JLabel("Select Tracks to Display");
            header.setBorder(BorderFactory.createEmptyBorder(12, 12, 0, 12));
            chooser.getContentPane().add(header, BorderLayout.NORTH);
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
        locator.setTrackNames(chosen);
        return !chosen.isEmpty();
    }

    /** The chooser's columns: the track's name, then its attributes in order of first appearance, at most 15. */
    static List<String> headings(FigWigFile file) {
        LinkedHashSet<String> names = new LinkedHashSet<>();
        names.add("name");
        for (int t = 0; t < file.nTracks(); t++) {
            names.addAll(file.trackAttributes(t).keySet());
        }
        List<String> headings = new ArrayList<>(names);
        return headings.size() > 15 ? headings.subList(0, 15) : headings;
    }

    static List<FileRecord> records(FigWigFile file) {
        List<FileRecord> records = new ArrayList<>(file.nTracks());
        for (int t = 0; t < file.nTracks(); t++) {
            Map<String, String> attributes = new LinkedHashMap<>();
            attributes.put("name", file.trackName(t));
            attributes.putAll(file.trackAttributes(t));
            records.add(new FileRecord(file.trackName(t), attributes));
        }
        return records;
    }
}
