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

package org.broad.igv.ui.action;

import org.broad.igv.encode.TrackChooser;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.encode.EncodeTrackChooserFactory;
import org.broad.igv.encode.FileRecord;

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
        OTHER
    }

    private static Logger log = LogManager.getLogger(BrowseEncodeAction.class);

    private static Map<String, Color> colors;

    static {
        colors = new HashMap<String, Color>();
        colors.put("H3K27AC", new Color(200, 0, 0));
        colors.put("H3K27ME3", new Color(200, 0, 0));
        colors.put("H3K36ME3", new Color(0, 0, 150));
        colors.put("H3K4ME1", new Color(0, 150, 0));
        colors.put("H3K4ME2", new Color(0, 150, 0));
        colors.put("H3K4ME3", new Color(0, 150, 0));
        colors.put("H3K9AC", new Color(100, 0, 0));
        colors.put("H3K9ME1", new Color(100, 0, 0));
    }

    static Set<String> sampleInfoAttributes = new HashSet<>(Arrays.asList(
            "dataType", "cell", "antibody", "lab", "Biosample", "AssayType", "Target"));

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
        String antibody = attributes.containsKey("antibody") ? attributes.get("antibody") : attributes.get("Target");
        if (antibody != null) {
            rl.setColor(colors.get(antibody.toUpperCase()));
        }
        for (Map.Entry<String, String> entry : attributes.entrySet()) {
            String value = entry.getValue();
            if (value != null && value.length() > 0 && sampleInfoAttributes.contains(entry.getKey())) {
                AttributeManager.getInstance().addAttribute(rl.getName(), entry.getKey(), value);
            }
        }
        rl.setMetadata(attributes);
        return rl;
    }

    /**
     * Return a friendly name for the track.  Unfortunately it is neccessary to hardcode certain attributes.
     *
     * @return
     */
    public String getTrackName(FileRecord record) {

        Map<String, String> attributes = record.getAttributes();
        StringBuffer sb = new StringBuffer();
        if (attributes.containsKey("cell")) sb.append(attributes.get("cell") + " ");
        if (attributes.containsKey("antibody")) sb.append(attributes.get("antibody") + " ");
        if (attributes.containsKey("dataType")) sb.append(attributes.get("dataType") + " ");
        if (attributes.containsKey("view")) sb.append(attributes.get("view") + " ");
        if (attributes.containsKey("replicate")) sb.append("rep " + attributes.get("replicate"));

        String trackName = sb.toString().trim();
        if (sb.length() == 0) {
            trackName = (new File(record.getPath())).getName();
        }

        return trackName;

    }
}
