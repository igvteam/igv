/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2017 Broad Institute
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
package org.broad.igv.ui.javafx.panel;

import javafx.beans.property.DoubleProperty;
import javafx.scene.layout.HBox;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.javafx.JavaFXUIUtilities;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.*;

// Intended as the rough equivalent of the DataPanelContainer class of the Swing UI.  Work in progress.
// Note: Not dealing with DnD yet.
public class DataPaneContainer extends HBox {
    private TrackRow trackRow = null;
    private List<DataPane> dataPanes = new ArrayList<DataPane>();


    public DataPaneContainer(TrackRow trackRow) {
        this.trackRow = trackRow;
        createDataPanes();
    }

    public void createDataPanes() {
        getChildren().clear();
        dataPanes.clear();

        for (ReferenceFrame f : FrameManager.getFrames()) {
            if (f.isVisible()) {
                DataPane dp = new DataPane(f, this);
                dp.backgroundProperty().bind(backgroundProperty());
                JavaFXUIUtilities.bindHeightToContainer(this, dp);
                getChildren().add(dp);
            }
        }
    }

    public DoubleProperty frameSpacingProperty() {
        return spacingProperty();
    }

    public Collection<TrackGroup> getTrackGroups() {
        return trackRow.getGroups();
    }

    // *** The following methods below this point copied over from TrackPanel as the functionality is the same. ***

    private void autoscale() {


        final Collection<Track> trackList = IGV.getInstance().getAllTracks();

        Map<String, List<Track>> autoscaleGroups = new HashMap<String, List<Track>>();

        for (Track track : trackList) {

            if (!track.isVisible()) continue;

            String asGroup = track.getAttributeValue(AttributeManager.GROUP_AUTOSCALE);
            if (asGroup != null) {
                if (!autoscaleGroups.containsKey(asGroup)) {
                    autoscaleGroups.put(asGroup, new ArrayList<Track>());
                }

                if (track instanceof MergedTracks) {
                    for (Track mt : ((MergedTracks) track).getMemberTracks()) {
                        // TODO: Is this a bug? Copied the code over like this, but seems it should be .add(mt)
                        autoscaleGroups.get(asGroup).add(track);
                    }
                } else {
                    autoscaleGroups.get(asGroup).add(track);
                }
            } else if (track.getAutoScale()) {

                if (track instanceof MergedTracks) {
                    for (Track mt : ((MergedTracks) track).getMemberTracks()) {
                        autoscaleGroup(Arrays.asList(mt));
                    }
                } else {
                    autoscaleGroup(Arrays.asList(track));
                }
            }

        }

        if (autoscaleGroups.size() > 0) {
            for (List<Track> tracks : autoscaleGroups.values()) {
                autoscaleGroup(tracks);
            }
        }
    }

    private void autoscaleGroup(List<Track> trackList) {


        List<ReferenceFrame> frames =
                FrameManager.isGeneListMode() ? FrameManager.getFrames() :
                        Arrays.asList(FrameManager.getDefaultFrame());


        List<Range> inViewRanges = new ArrayList<Range>();

        synchronized (trackList) {
            for (Track track : trackList) {
                if (track instanceof ScalableTrack) {
                    for (ReferenceFrame frame : frames) {
                        // TODO: port to JavaFX
//                        Range range = ((ScalableTrack) track).getInViewRange(frame);
//                        if (range != null) {
//                            inViewRanges.add(range);
//                        }
                    }
                }
            }

            if (inViewRanges.size() > 0) {

                Range inter = computeScale(inViewRanges);

                for (Track track : trackList) {

                    DataRange dr = track.getDataRange();
                    float min = Math.min(0, inter.min);
                    float base = Math.max(min, dr.getBaseline());
                    float max = inter.max;
                    // Pathological case where min ~= max  (no data in view)
                    if (max - min <= (2 * Float.MIN_VALUE)) {
                        max = min + 1;
                    }

                    DataRange newDR = new DataRange(min, base, max, dr.isDrawBaseline());
                    newDR.setType(dr.getType());
                    track.setDataRange(newDR);

                }
            }
        }
    }

    public static Range computeScale(List<Range> ranges) {

        float min = 0;
        float max = 0;

        if (ranges.size() > 0) {
            max = ranges.get(0).max;
            min = ranges.get(0).min;

            for (int i = 1; i < ranges.size(); i++) {

                Range r = ranges.get(i);
                max = Math.max(r.max, max);
                min = Math.min(r.min, min);

            }
        }

        return new Range(min, max);
    }
}
