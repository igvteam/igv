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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

//~--- JDK imports ------------------------------------------------------------

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.List;
import java.util.*;

/**
 * Container for a group of tracks.  Behaves as a single unit when sorting
 * by region score.
 *
 * @author jrobinso
 */
public class TrackPanelHelper {

    private static Logger log = LogManager.getLogger(TrackPanelHelper.class);

    public static List<Track> sortByAttributes(
            List<Track> tracks,
            final String[] attributeNames,
            final boolean[] ascending) {


        if ((tracks != null) && !tracks.isEmpty()) {
            List<Track> allTracks = new ArrayList<Track>(tracks);
            try {
                Comparator<Track> comparator = new TrackAttributeComparator(attributeNames, ascending);

                // Step 1, remove non-sortable tracks and remember position
                List<Track> nonsortableTracks = new ArrayList<Track>();
                Map<Track, Integer> trackIndices = new HashMap<Track, Integer>();
                for (int i = tracks.size() - 1; i >= 0; i--) {
                    if (!tracks.get(i).isSortable()) {
                        Track t = tracks.remove(i);
                        nonsortableTracks.add(t);
                        trackIndices.put(t, i);

                    }
                }

                // Step 2, sort "sortable" tracks
                Collections.sort(tracks, comparator);

                // Step 2.5, internal sort by sample attributes for variant tracks.  This is ugly but neccessary as
                // variant tracks are implemented as monoliths, with sample rows internal to the track.
                for (Track t : allTracks) {
                    if (t instanceof org.broad.igv.variant.VariantTrack) {
                        ((org.broad.igv.variant.VariantTrack) t).sortSamples(new SampleAttributeComparator(attributeNames, ascending));
                    }
                }

                // Step 3, put non-sortable tracks back in original order
                if (nonsortableTracks.size() > 0) {
                    for (int i = nonsortableTracks.size() - 1; i >= 0; i--) {
                        Track t = nonsortableTracks.get(i);
                        int index = trackIndices.get(t);
                        tracks.add(index, t);
                    }
                }
            } catch (Exception e) {
                log.error("Error sorting tracks by attribute", e);
                tracks = allTracks;
            }
        }
        return tracks;
    }


    /**
     * Sorts the entire group, including tracks for the given score type as well as other tracks, by
     * the specified sample order
     */
    public static List<Track> sortByRegionScore(
            List<Track> tracks,
            final RegionScoreType type,
            List<String> sortedSamples) {

        // Step 1,  remove non-sortable tracks and remember position
        List<Track> unsortableTracks = new ArrayList();
        Map<Track, Integer> trackIndeces = new HashMap();
        for (int i = tracks.size() - 1; i >= 0; i--) {
            if (!tracks.get(i).isSortable()) {
                Track t = tracks.remove(i);
                unsortableTracks.add(t);
                trackIndeces.put(t, i);

            }
        }

        List<Track> tracksWithScore = new ArrayList();
        List<Track> otherTracks = new ArrayList();
        for (Track t : tracks) {
            if (t.isRegionScoreType(type)) {
                tracksWithScore.add(t);
            } else {
                otherTracks.add(t);
            }
        }

        sortBySampleOrder(tracksWithScore, sortedSamples);
        sortBySampleOrder(otherTracks, sortedSamples);

        tracks.clear();
        tracks.addAll(tracksWithScore);
        tracks.addAll(otherTracks);

        // Step 3, put unortable tracks back in original order
        if (unsortableTracks.size() > 0) {
            for (int i = unsortableTracks.size() - 1; i >= 0; i--) {
                Track t = unsortableTracks.get(i);
                int index = trackIndeces.get(t);
                tracks.add(index, t);
            }
        }
        return tracks;
    }

    /**
     * @param sortedSamples
     */
    private static List<Track> sortBySampleOrder(List<Track> tracks,
                                                 List<String> sortedSamples) {
        if ((tracks != null) && (sortedSamples != null) && !tracks.isEmpty()) {

            // Create a rank hash.  Loop backwards so that the lowest index for an attribute
            final HashMap<String, Integer> rankMap = new HashMap(sortedSamples.size() * 2);
            for (int i = sortedSamples.size() - 1; i >= 0; i--) {
                rankMap.put(sortedSamples.get(i), i);
            }
            // Comparator for sorting in ascending order
            Comparator<Track> c = new Comparator<Track>() {

                public int compare(Track t1, Track t2) {
                    String a1 = t1.getSample(); //t1.getAttributeValue(attributeId);
                    String a2 = t2.getSample(); //t2.getAttributeValue(attributeId);
                    Integer r1 = ((a1 == null) ? null : rankMap.get(a1));
                    Integer r2 = ((a2 == null) ? null : rankMap.get(a2));
                    if ((r1 == null) && (r2 == null)) {
                        return 0;
                    } else if (r1 == null) {
                        return 1;
                    } else if (r2 == null) {
                        return -1;
                    } else {
                        return r1.intValue() - r2.intValue();
                    }

                }
            };

            Collections.sort(tracks, c);
        }
        return tracks;
    }

    /**
     * @param trackIds
     */
    public static List<Track> sortByList(List<Track> tracks, List<String> trackIds) {

        final Map<String, Integer> trackPositions = new HashMap();
        for (int i = 0; i < trackIds.size(); i++) {
            trackPositions.put(trackIds.get(i), i);
        }
        Comparator c = new Comparator<Track>() {
            public int compare(Track t1, Track t2) {
                String id1 = t1.getId();
                int p1 = trackPositions.containsKey(id1) ? trackPositions.get(id1) : Integer.MAX_VALUE;
                String id2 = t2.getId();
                int p2 = trackPositions.containsKey(id2) ? trackPositions.get(id2) : Integer.MAX_VALUE;
                return p1 - p2;
            }
        };
        Collections.sort(tracks, c);

        return tracks;
    }


    private abstract static class AttributeComparator<T> implements Comparator<T> {

        private final String[] attributeNames;
        private final boolean[] ascending;

        AttributeComparator(String[] attributeNames, boolean[] ascending) {
            assert attributeNames.length == ascending.length;
            this.attributeNames = attributeNames;
            this.ascending = ascending;
        }

        protected abstract String getAttributeValue(T track, String attName);

        public int compare(T t1, T t2) {
            // Loop through the attributes in order (primary, secondary, tertiary, ...).  The
            // first attribute to yield a non-zero comparison wins
            for (int i = 0; i < attributeNames.length; i++) {
                String attName = attributeNames[i];

                if (attName != null) {
                    String value1 = getAttributeValue(t1, attName);
                    String value2 = getAttributeValue(t2, attName);

                    boolean isNumeric = AttributeManager.getInstance().isNumeric(attName);

                    int c = 0;
                    if (isNumeric) {
                        double d1;
                        try {
                            d1 = Double.parseDouble(value1);
                        } catch (NumberFormatException e) {
                            d1 = Double.MIN_VALUE;
                        }
                        double d2;
                        try {
                            d2 = Double.parseDouble(value2);
                        } catch (NumberFormatException e) {
                            d2 = Double.MIN_VALUE;
                        }
                        c = Double.compare(d1, d2);
                    } else {
                        c = value1.compareTo(value2);
                    }

                    if (c != 0) {
                        return ascending[i] ? c : -c;
                    }

                }
            }

            // All compares are equal
            return 0;
        }
    }

    private static class TrackAttributeComparator extends AttributeComparator<Track> {

        public TrackAttributeComparator(String[] attributeNames, boolean[] ascending) {
            super(attributeNames, ascending);
        }

        protected String getAttributeValue(Track track, String attName) {
            String value = track.getAttributeValue(attName);

            if (value == null) {
                value = "";
            }

            return value.toLowerCase();
        }
    }


    /**
     * Sort samples by attribute value; logic copied wholesale from
     * AttributeComparator, probably some refactoring could be done
     * to minimize the duplicated code
     */
    private static class SampleAttributeComparator extends AttributeComparator<String> {

        public SampleAttributeComparator(String[] attributeNames, boolean[] ascending) {
            super(attributeNames, ascending);
        }

        protected String getAttributeValue(String sample, String attName) {
            String value = AttributeManager.getInstance().getAttribute(sample, attName);

            if (value == null) {
                value = "";
            }

            return value.toLowerCase();
        }

    }


}