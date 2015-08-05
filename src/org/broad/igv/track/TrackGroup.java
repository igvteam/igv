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
package org.broad.igv.track;

//~--- JDK imports ------------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * Container for a group of tracks.  Behaves as a single unit when sorting
 * by region score.
 *
 * @author jrobinso
 */
public class TrackGroup {

    private static Logger log = Logger.getLogger(TrackGroup.class);

    /**
     * Key used to group tracks (e.g. SAMPLE_ID).
     */
    private String name;

    private boolean drawBorder = true;

    private List<Track> tracks;

    private boolean selected;


    public TrackGroup() {
        this("");
    }


    public TrackGroup(String name) {
        this.name = name;
        tracks = Collections.synchronizedList(new ArrayList<Track>());
    }


    public boolean contains(Track track) {
        return tracks.contains(track);
    }

    public List<Track> getTracks() {
        return tracks;
    }


    public int indexOf(Track track) {
        return tracks.indexOf(track);
    }


    public int size() {
        return tracks.size();
    }


    public void add(Track track) {
        if (track == null) {
            log.info("Attempt to add null track");
        } else {
            log.debug("Adding track " + track.getName());
            tracks.add(track);
        }
    }

    public void add(int pos, Track track) {
        if (track == null) {
            log.info("Attempt to add null track");
        } else {
            tracks.add(pos, track);
        }
    }

    public void addAll(Collection<Track> trackList) {
        tracks.addAll(trackList);
    }

    public void addAll(int index, Collection<Track> trackList) {
        tracks.addAll(index, trackList);
    }

    public void remove(Track track) {
        tracks.remove(track);
    }

    /**
     * Return a composite score for the entire group.  For now use the maximum track
     * score.   Note that scores for tracks not appropriate to the score type will
     * return -Float.MAX, so they are effectively ignored.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @param frameName
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        float score = -Float.MAX_VALUE;
        for (Track track : tracks) {
            if (track.isVisible()) {
                score = Math.max(score, track.getRegionScore(chr, start, end, zoom, type, frameName));

            }
        }
        return score;
    }


    public String getName() {
        return name;
    }

    public boolean isDrawBorder() {
        return drawBorder;
    }


    public void setDrawBorder(boolean drawBorder) {
        this.drawBorder = drawBorder;
    }

    public boolean isVisible() {
        for (Track t : tracks) {
            if ((t != null) && t.isVisible()) {
                return true;
            }
        }
        return false;
    }

    public int getHeight() {

        int h = 0;
        for (Track track : tracks) {
            if (track != null && track.isVisible()) {
                h += track.getHeight();
            }
        }
        return h;
    }

    public void renderName(Graphics2D g2D, Rectangle rect, boolean isSelected) {

        // Calculate fontsize
        int fontSize = PreferenceManager.getInstance().getAsInt(PreferenceManager.DEFAULT_FONT_SIZE);

        Font font = FontManager.getFont(Font.BOLD, fontSize);
        g2D.setFont(font);

        GraphicUtils.drawWrappedText(getName(), rect, g2D, true);

    }

    /**
     * Sort tracks by the array of attribute names.
     *
     * @param attributeNames
     * @param ascending
     */

    public void sortByAttributes(final String[] attributeNames,
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

    }


    /**
     * Sorts the entire group, including tracks for the given score type as well as other tracks, by
     * the specified sample order
     */
    public void sortGroup(final RegionScoreType type,
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

        List<Track> tracksWithScore = new ArrayList(getTracks().size());
        List<Track> otherTracks = new ArrayList(getTracks().size());
        for (Track t : getTracks()) {
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

    }

    /**
     * @param sortedSamples
     */
    private void sortBySampleOrder(List<Track> tracks,
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

    }

    /**
     * @param trackIds
     */
    public void sortByList(List<String> trackIds) {

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
    }

    public void removeTracks(Collection<? extends Track> tracksToRemove) {
        tracks.removeAll(tracksToRemove);
    }


    /**
     * If this group contains the targetTrack, insert the selectedTracks collection either before or after
     * the target and return true.   Otherwise return false.
     *
     * @param selectedTracks
     * @param targetTrack
     * @param before
     */
    public boolean moveSelectedTracksTo(Collection<? extends Track> selectedTracks,
                                        Track targetTrack,
                                        boolean before) {

        int index = (targetTrack == null ? tracks.size() : tracks.indexOf(targetTrack));
        if (index < 0) {
            return false;
        }

        if (!before) {
            index = index + 1;
        }

        // 1. Divdide the target list up into 2 parts, one before the index and one after
        List<Track> beforeList = new ArrayList(tracks.subList(0, index));
        List<Track> afterList = new ArrayList(tracks.subList(index, tracks.size()));

        // 2.  Remove the selected tracks from anywhere they occur
        beforeList.removeAll(selectedTracks);
        afterList.removeAll(selectedTracks);

        // 3. Now insert the selected tracks
        tracks.clear();
        tracks.addAll(beforeList);
        tracks.addAll(selectedTracks);
        tracks.addAll(afterList);

        return true;

    }

    public boolean isSelected() {
        return selected;
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
        IGV igv = IGV.getInstance();
        igv.clearSelections();
        igv.setTrackSelections(new HashSet(tracks));

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