package org.broad.igv.ui;

import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.*;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.*;

public class Autoscaler {

    public static void autoscale(Collection<Track> trackList) {
        Map<String, List<Track>> autoscaleGroups = new HashMap<>();
        for (Track track : trackList) {
            if (!track.isVisible()) continue;

            String asGroup = track.getAttributeValue(AttributeManager.GROUP_AUTOSCALE);
            if (asGroup != null) {
                if (!autoscaleGroups.containsKey(asGroup)) {
                    autoscaleGroups.put(asGroup, new ArrayList<Track>());
                }
                if (track instanceof MergedTracks) {
                    for (Track mt : ((MergedTracks) track).getMemberTracks()) {
                        autoscaleGroups.get(asGroup).add(mt);
                    }
                } else {
                    autoscaleGroups.get(asGroup).add(track);
                }
            } else if (track.getAutoScale()) {
                if (track instanceof MergedTracks) {
                    List<Track> memberTracks = new ArrayList(((MergedTracks) track).getMemberTracks());
                    autoscaleGroup(memberTracks);

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

    private static void autoscaleGroup(List<Track> trackList) {
        List<ReferenceFrame> frames =
                FrameManager.isGeneListMode() ? FrameManager.getFrames() :
                        Arrays.asList(FrameManager.getDefaultFrame());

        List<Range> inViewRanges = Collections.synchronizedList(new ArrayList<>());
        for (Track track : trackList) {
            if (track instanceof ScalableTrack) {
                for (ReferenceFrame frame : frames) {
                    Range r = ((ScalableTrack) track).getInViewRange(frame);
                    if (r != null) {
                        inViewRanges.add(r);
                    }
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
