/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.session;


import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.lists.GeneList;
import org.broad.igv.track.TrackManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ObserverForObject;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author eflakes
 */
public class Session {
    private static Logger log = Logger.getLogger(Session.class);

    public static String[] PREF_KEYS = {
            PreferenceManager.DISPLAY_OVERLAY_TRACKS_KEY,
            PreferenceManager.COLOR_OVERLAY_KEY,
            PreferenceManager.OVERLAY_ATTRIBUTE_KEY

    };

    private static int version = 4;

    private TrackManager trackManager;

    private String filePath;
    private String groupTracksBy;
    private ReferenceFrame referenceFrame = FrameManager.getDefaultFrame();
    private TrackFilter filter;
    private HashMap<String, String> preferences;
    private HashMap<TrackType, ContinuousColorScale> colorScales;

    double [] dividerFractions = null;

    private History history;

    /**
     * Map of chromosome -> regions of interest
     */
    private Map<String, Collection<RegionOfInterest>> regionsOfInterest;
    //An Observable that notifies observers of changes to the regions of interest.  Its
    //setChangedAndNotify() method should be called after any regions change.
    private ObserverForObject<Map<String, Collection<RegionOfInterest>>> regionsOfInterestObservable;

    private GeneList currentGeneList;


    public Session(String filePath) {
        log.debug("New session");

        this.filePath = filePath;



        regionsOfInterest = new LinkedHashMap<String, Collection<RegionOfInterest>>();
        regionsOfInterestObservable =
                new ObserverForObject<Map<String, Collection<RegionOfInterest>>>(regionsOfInterest);

        preferences = new HashMap<String, String>();
        colorScales = new HashMap<TrackType, ContinuousColorScale>();
        history = new History(100);

        boolean resetRequired = FrameManager.getFrames().size() > 1;
        setCurrentGeneList(null);
        if (resetRequired) {
            IGV.getInstance().resetFrames();
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

    public String getSessionVersion() {
        return String.valueOf(version);
    }

    /**
     * @return the absolute path to the file associated with this session
     */
    public String getPath() {
        return filePath;
    }

    public void setPreference(String key, String value) {
        preferences.put(key, value);
    }


    public String getPreference(String key) {
        return preferences.get(key);
    }


    public void setColorScale(TrackType trackType, ContinuousColorScale colorScale) {
        colorScales.put(trackType, colorScale);
    }

    public ContinuousColorScale getColorScale(TrackType trackType) {
        if (colorScales.containsKey(trackType)) {
            return colorScales.get(trackType);
        } else {
            return PreferenceManager.getInstance().getColorScale(trackType);
        }
    }

    public boolean getDisplayOverlayTracks() {
        final String key = PreferenceManager.DISPLAY_OVERLAY_TRACKS_KEY;
        if (preferences.containsKey(key)) {
            try {
                return Boolean.parseBoolean(preferences.get(key));
            }
            catch (Exception e) {
                log.error("Error converting boolean preference " + key + "=" + preferences.get(key));
            }
        }
        return PreferenceManager.getInstance().getAsBoolean(PreferenceManager.DISPLAY_OVERLAY_TRACKS_KEY);

    }

    public boolean getColorOverlay() {
        final String key = PreferenceManager.COLOR_OVERLAY_KEY;
        if (preferences.containsKey(key)) {
            try {
                return Boolean.parseBoolean(preferences.get(key));
            }
            catch (Exception e) {
                log.error("Error converting boolean preference " + key + "=" + preferences.get(key));
            }
        }
        return PreferenceManager.getInstance().getAsBoolean(PreferenceManager.COLOR_OVERLAY_KEY);

    }

    public String getOverlayAttribute() {
        final String key = PreferenceManager.OVERLAY_ATTRIBUTE_KEY;
        if (preferences.containsKey(key)) {
            return preferences.get(key);

        }
        return PreferenceManager.getInstance().get(PreferenceManager.OVERLAY_ATTRIBUTE_KEY);
    }

    public String getTrackAttributeName() {
        final String key = PreferenceManager.TRACK_ATTRIBUTE_NAME_KEY;
        if (preferences.containsKey(key)) {
            return preferences.get(key);

        }
        return PreferenceManager.getInstance().get(PreferenceManager.TRACK_ATTRIBUTE_NAME_KEY);
    }

    public String getLocusString() {
        if(getReferenceFrame().getChrName().equals(Globals.CHR_ALL)) {
            return Globals.CHR_ALL;
        }
        ReferenceFrame.Range range = getReferenceFrame().getCurrentRange();
        String startStr = String.valueOf(range.getStart());
        String endStr = String.valueOf(range.getEnd());
        String position = range.getChr() + ":" + startStr + "-" + endStr;
        return position;
    }

    public void setLocus(String locusString) {
        try {
            Locus locus = getLocus(locusString);
            if (locus == null) {
                Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
                if (genome != null) {
                    referenceFrame.setChrName(genome.getHomeChromosome());
                }
            } else {
               referenceFrame.jumpTo(locus.chr, locus.start, locus.end);
            }
            getHistory().push(locusString, referenceFrame.getZoom());
        }
        catch (Exception e) {
            MessageUtils.showMessage("Error setting locus string: " + locusString + " (" + e.toString() + ")");
            log.error("Error setting locus string", e);
        }
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
     * Removes the regions of interest from the current chromosome.  returns global success/failure.
     * This method to remove multiple regions at once exists because, if you do them one at a time,
     * it throws the RegionNavigatorDialog table rows off.
     *
     * @param rois
     * @return whether all specified regions were found for removal
     */
    public boolean removeRegionsOfInterest(Collection<RegionOfInterest> rois) {
        String chr = FrameManager.getDefaultFrame().getChrName();
        Collection<RegionOfInterest> roiList = regionsOfInterest.get(chr);
        boolean result = true;
        if (roiList != null) {
            for (RegionOfInterest roi : rois) {
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

    public void setFilePath(String filePath) {
        this.filePath = filePath;
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
        this.currentGeneList = currentGeneList;
        FrameManager.resetFrames(currentGeneList);
    }

    /**
     * Add a gene to the current list.  This is a transient change (not saved to disk)
     *
     * @param gene
     */
    public void addGene(String gene) {
        if (getCurrentGeneList() != null) {
            getCurrentGeneList().add(gene);
            setCurrentGeneList(getCurrentGeneList());
        }
    }

    public int getVersion() {
        return version;
    }

    public void setVersion(int version) {
        this.version = version;
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


    static Locus getLocus(String searchString) {

        if (searchString != null) {
            String chr;
            int[] startEnd;
            int colon = searchString.indexOf(":");

            if (colon > 0) {

                // The chromosome is that portion of the search string up to the colon.
                chr = searchString.substring(0, colon);
                String posString = searchString.substring(colon).replace(":", "");
                startEnd = getStartEnd(posString);

                if (startEnd != null) {
                    return new Locus(chr, startEnd[0], startEnd[1]);
                }
            } else {

                // No chromosome delimiter (color),  The search string is either chromosome name
                // or a locus in the current chromosome.
                if (searchString.contains("-")) {

                    // Presense of a dash indicates this is a locus string in the current chromosome
                    startEnd = getStartEnd(searchString);
                    if (startEnd != null) {
                        return new Locus(null, startEnd[0], startEnd[1]);
                    }
                } else {
                    Feature feature = FeatureDB.getFeature(searchString.toUpperCase().trim());
                    if (feature != null) {
                        return new Locus(feature.getChr(), feature.getStart(), feature.getEnd());
                    } else {

                        Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
                        if (genome.getChromosome(searchString) != null) {
                            // No dash, this is either a chromosome or an unkown search string
                            return new Locus(searchString, -1, -1);
                        } else {
                            return null;
                        }
                    }
                }
            }
        }
        return null;
    }


    /**
     * Return the start and end positions as a 2 element array for the input
     * position string.  UCSC conventions  are followed for coordinates,
     * specifically the internal representation is "zero" based (first base is
     * numbered 0) but the display representation is "one" based (first base is
     * numbered 1).   Consequently 1 is substracted from the parsed positions
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
    public ObserverForObject<Map<String, Collection<RegionOfInterest>>> getRegionsOfInterestObservable() {
        return regionsOfInterestObservable;
    }
}


