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

package org.broad.igv.session;

import org.broad.igv.Globals;
import org.broad.igv.event.IGVEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.event.ViewChange;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.lists.GeneList;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.sam.InsertionManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ObservableForObject;

import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author eflakes
 */
public class Session implements IGVEventObserver {

    private static Logger log = LogManager.getLogger(Session.class);

    //This doesn't mean genelist or not, the same way it does in FrameManager
    public enum GeneListMode {
        NORMAL, CURSOR
    }

    private int version;
    private String path;
    private String groupTracksBy;
    public boolean expandInsertions = false; //false;
    private int nextAutoscaleGroup;
    private ReferenceFrame referenceFrame = FrameManager.getDefaultFrame();
    private TrackFilter filter;
    private HashMap<String, String> preferences;
    private HashMap<TrackType, ContinuousColorScale> trackColorScales;
    private boolean removeEmptyPanels = false;
    double[] dividerFractions = null;

    /**
     * Attribute used to group tracks.  Normally "null".  Set from the "Tracks" menu.
     */
    private String groupByAttribute = null;



    private History history;

    /**
     * Map of chromosome -> regions of interest
     */
    private Map<String, Collection<RegionOfInterest>> regionsOfInterest;
    //An Observable that notifies observers of changes to the regions of interest.  Its
    //setChangedAndNotify() method should be called after any regions change.
    private ObservableForObject<Map<String, Collection<RegionOfInterest>>> regionsOfInterestObservable;

    private GeneList currentGeneList;
    private GeneListMode geneListMode = GeneListMode.NORMAL;
    private Set<String> hiddenAttributes;
    private String locus;

    public Session(String path) {
        init(path);
    }

    private void init(String path) {
        this.path = path;
        this.nextAutoscaleGroup = 1;
        this.groupByAttribute = null;
        this.regionsOfInterest = new LinkedHashMap<>();
        this.regionsOfInterestObservable = new ObservableForObject<>(regionsOfInterest);
        this.preferences = new HashMap<>();
        this.trackColorScales = new HashMap<>();
        this.hiddenAttributes = null;
        this.history = new History(100);
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
    }

    public void reset(String path) {
        init(path);
        setCurrentGeneList(null);
        if (FrameManager.getFrames().size() > 1) {
            IGV.getInstance().resetFrames();
        }
        for(ReferenceFrame frame : FrameManager.getFrames()) {
            frame.setExpandedInsertion(null);
        }
        InsertionManager.getInstance().clear();

    }


    public void receiveEvent(IGVEvent event) {
        if (event instanceof ViewChange) {
            ViewChange e = (ViewChange) event;
            if (e.recordHistory()) {
                recordHistory();
            }
        } else {
            log.warn("Unknown event type: " + event.getClass());
        }
    }

    public void clearDividerLocations() {
        dividerFractions = null;
    }

    public void setDividerFractions(double[] divs) {
        this.dividerFractions = divs;
    }

    public double[] getDividerFractions() {
        return dividerFractions;
    }

    public void recordHistory() {
        final ReferenceFrame defaultFrame = FrameManager.getDefaultFrame();
        history.push(defaultFrame.getFormattedLocusString(), defaultFrame.getZoom());
    }

    public void clearHistory() {
        history.clear();
    }

    public String getSessionVersion() {
        return String.valueOf(version);
    }

    /**
     * @return the absolute path to the file associated with this session.  This can be null if session was not
     * initalized from a file.
     */
    public String getPath() {
        return path;
    }

    public String getGroupByAttribute() {
        return groupByAttribute;
    }

    public void setGroupByAttribute(String groupByAttribute) {
        this.groupByAttribute = groupByAttribute;
    }

    /**
     * Set a preference value for this session only
     *
     * @param key
     * @param value
     */
    public void setPreference(String key, String value) {
        preferences.put(key, value);
    }

    public void removePreference(String key) {
        preferences.remove(key);
    }

    public void setColorScale(TrackType trackType, ContinuousColorScale colorScale) {
        trackColorScales.put(trackType, colorScale);
    }

    public ContinuousColorScale getColorScale(TrackType trackType) {
        if (trackColorScales.containsKey(trackType)) {
            return trackColorScales.get(trackType);
        } else {
            return PreferencesManager.getPreferences().getColorScale(trackType);
        }
    }


    public boolean getPreferenceAsBoolean(String key) {
        if (preferences.containsKey(key)) {
            try {
                return Boolean.parseBoolean(preferences.get(key));
            } catch (Exception e) {
                log.error("Error converting boolean preference " + key + "=" + preferences.get(key));
            }
        }
        return PreferencesManager.getPreferences().getAsBoolean(key);

    }

    public String getPreference(String key) {
        if (preferences.containsKey(key)) {
            return (preferences.get(key));
        } else {
            return PreferencesManager.getPreferences().get(key);
        }
    }

    /**
     * @param key
     * @param def Default value. Not saved
     * @return Preference if found, or else default value
     * @see IGVPreferences#getPersistent(String, String)
     */
    public String getPersistent(String key, String def) {
        if (preferences.containsKey(key)) {
            return (preferences.get(key));
        } else {
            return PreferencesManager.getPreferences().getPersistent(key, def);
        }
    }

    public boolean getOverlayMutationTracks() {
        final String key = OVERLAY_MUTATION_TRACKS;
        if (preferences.containsKey(key)) {
            try {
                return Boolean.parseBoolean(preferences.get(key));
            } catch (Exception e) {
                log.error("Error converting boolean preference " + key + "=" + preferences.get(key));
            }
        }
        return PreferencesManager.getPreferences().getAsBoolean(OVERLAY_MUTATION_TRACKS);

    }

    public boolean getColorOverlay() {
        final String key = COLOR_MUTATIONS;
        if (preferences.containsKey(key)) {
            try {
                return Boolean.parseBoolean(preferences.get(key));
            } catch (Exception e) {
                log.error("Error converting boolean preference " + key + "=" + preferences.get(key));
            }
        }
        return PreferencesManager.getPreferences().getAsBoolean(COLOR_MUTATIONS);

    }

    public String getOverlayAttribute() {
        final String key = OVERLAY_ATTRIBUTE_KEY;
        if (preferences.containsKey(key)) {
            return preferences.get(key);

        }
        return PreferencesManager.getPreferences().get(OVERLAY_ATTRIBUTE_KEY);
    }

    public String getTrackAttributeName() {
        final String key = TRACK_ATTRIBUTE_NAME_KEY;
        if (preferences.containsKey(key)) {
            return preferences.get(key);
        }
        return PreferencesManager.getPreferences().get(TRACK_ATTRIBUTE_NAME_KEY);
    }


    public String getLocusString() {
        if (getReferenceFrame().getChrName().equals(Globals.CHR_ALL)) {
            return Globals.CHR_ALL;
        }
        Range range = getReferenceFrame().getCurrentRange();
        String startStr = String.valueOf(range.getStart());
        String endStr = String.valueOf(range.getEnd());
        String position = range.getChr() + ":" + startStr + "-" + endStr;
        return position;
    }

    public void setLocus(String locus) {
        this.locus = locus;
    }

    public String getLocus() {
        return locus;
    }

    public TrackFilter getFilter() {
        return filter;
    }

    public void setFilter(TrackFilter filter) {
        this.filter = filter;
    }

    public String getGroupTracksBy() {
        return groupTracksBy;
    }

    public void setGroupTracksBy(String groupTracksBy) {
        this.groupTracksBy = groupTracksBy;
    }

    public ReferenceFrame getReferenceFrame() {
        return referenceFrame;
    }

    /**
     * WARNING! This method should never be used for update of the regions collection.
     * If this gets abused, this method could be changed to create a new Collection so that updates
     * have no effect on the real one.
     *
     * @param chr
     * @return
     */
    public Collection<RegionOfInterest> getRegionsOfInterest(String chr) {
        if (chr.equals(Globals.CHR_ALL)) {
            return getAllRegionsOfInterest();
        } else {
            return regionsOfInterest.get(chr);
        }
    }


    public Collection<RegionOfInterest> getAllRegionsOfInterest() {
        ArrayList<RegionOfInterest> roiList = new ArrayList<RegionOfInterest>();
        for (Collection<RegionOfInterest> roi : regionsOfInterest.values()) {
            roiList.addAll(roi);
        }
        return roiList;
    }

    /**
     * Removes the regions of interest.  returns global success/failure.
     * This method to remove multiple regions at once exists because, if you do them one at a time,
     * it throws the RegionNavigatorDialog table rows off.
     *
     * @param rois
     * @return whether all specified regions were found for removal
     */
    public boolean removeRegionsOfInterest(Collection<RegionOfInterest> rois) {
        boolean result = true;

        for (RegionOfInterest roi : rois) {
            Collection<RegionOfInterest> roiList = regionsOfInterest.get(roi.getChr());
            if (roiList != null) {
                result = result && roiList.remove(roi);
            }
        }

        //notify all observers that regions have changed.
        regionsOfInterestObservable.setChangedAndNotify();
        return result;
    }


    public void addRegionOfInterestWithNoListeners(RegionOfInterest roi) {
        String chr = roi.getChr();
        Collection<RegionOfInterest> roiList = regionsOfInterest.get(chr);
        if (roiList == null) {
            roiList = new ArrayList<RegionOfInterest>();
            regionsOfInterest.put(chr, roiList);
        }
        roiList.add(roi);

        //notify all observers that regions have changed.
        regionsOfInterestObservable.setChangedAndNotify();
    }

    public void clearRegionsOfInterest() {
        if (regionsOfInterest != null) {
            regionsOfInterest.clear();
        }
        //notify all observers that regions have changed.
        regionsOfInterestObservable.setChangedAndNotify();
    }

    public void setPath(String path) {
        this.path = path;
    }

    public History getHistory() {
        return history;
    }

    public List<History.Entry> getAllHistory() {
        return getHistory().getAllHistory();
    }

    public GeneList getCurrentGeneList() {
        return currentGeneList;
    }

    public void setCurrentGeneList(GeneList currentGeneList) {

        boolean frameReset = (currentGeneList != null || FrameManager.isGeneListMode());
        this.currentGeneList = currentGeneList;

        if (frameReset) {
            FrameManager.resetFrames(currentGeneList);
        }
    }

    public int getVersion() {
        return version;
    }

    public void setVersion(int version) {
        this.version = version;
    }


    public void setNextAutoscaleGroup(int nextAutoscaleGroup) {
        this.nextAutoscaleGroup = nextAutoscaleGroup;
    }

    public int getNextAutoscaleGroup() {
        return this.nextAutoscaleGroup++;
    }

    /**
     * Return a set containing names of attributes explicitly marked hidden.  If no attributes have been explicitly
     * marked return the default set.
     *
     * @return
     */
    public Set<String> getHiddenAttributes() {
        if (hiddenAttributes == null) {
            return (PreferencesManager.getPreferences().getAsBoolean(SHOW_DEFAULT_TRACK_ATTRIBUTES))  ?
                    Collections.emptySet() :
                    new HashSet<>(AttributeManager.defaultTrackAttributes);
        } else {
            return hiddenAttributes;
        }
    }

    public void setHiddenAttributes(Set<String> attributes) {
        this.hiddenAttributes = attributes;

    }

    public boolean isRemoveEmptyPanels() {
        return removeEmptyPanels;
    }

    public void setRemoveEmptyPanels(boolean removeEmptyPanels) {
        this.removeEmptyPanels = removeEmptyPanels;
    }


    static class Locus {
        String chr;
        int start;
        int end;

        Locus(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }

    }


    /**
     * Return the start and end positions as a 2 element array for the input
     * position string.  UCSC conventions  are followed for coordinates,
     * specifically the internal representation is "zero" based (first base is
     * numbered 0) but the display representation is "one" based (first base is
     * numbered 1).   Consequently, 1 is subtracted from the parsed positions
     */
    private static int[] getStartEnd(String posString) {
        try {
            String[] posTokens = posString.split("-");
            String startString = posTokens[0].replaceAll(",", "");
            int start = Math.max(0, Integer.parseInt(startString)) - 1;

            // Default value for end

            int end = start + 1;
            if (posTokens.length > 1) {
                String endString = posTokens[1].replaceAll(",", "");
                end = Integer.parseInt(endString);
            }

            if (posTokens.length == 1 || (end - start) < 10) {
                int center = (start + end) / 2;
                start = center - 20;
                end = center + 20;
            } else {
                String endString = posTokens[1].replaceAll(",", "");

                // Add 1 bp to end position t make it "inclusive"
                end = Integer.parseInt(endString);
            }

            return new int[]{Math.min(start, end), Math.max(start, end)};
        } catch (NumberFormatException numberFormatException) {
            return null;
        }

    }

    /**
     * Allows access to the Observable that notifies of changes to the regions of interest
     *
     * @return
     */
    public ObservableForObject<Map<String, Collection<RegionOfInterest>>> getRegionsOfInterestObservable() {
        return regionsOfInterestObservable;
    }
}


