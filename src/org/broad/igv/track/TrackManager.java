/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDescriptor;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.ui.MessageCollection;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.vcf.VCFTrack;
import org.broad.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.List;
import java.util.zip.GZIPInputStream;


/**
 * A helper class to manage tracks and groups.
 */

public class TrackManager {

    private static Logger log = Logger.getLogger(TrackManager.class);

    /**
     * Owning object
     */
    IGV igv;

    private TrackLoader loader;

    private final Map<String, TrackPanelScrollPane> trackPanelScrollPanes = new Hashtable();

    /**
     * Attribute used to group tracks.  Normally "null".  Set from the "Tracks" menu.
     */
    private String groupByAttribute = null;

    /**
     * The gene track for the current genome, rendered in the FeaturePanel
     */
    private Track geneTrack;

    /**
     * The sequence track for the current genome
     */
    private SequenceTrack sequenceTrack;

    private Map<String, List<Track>> overlayTracksMap = new HashMap();

    private Set<TrackType> loadedTypes = new HashSet();
    public static final String DATA_PANEL_NAME = "DataPanel";
    public static final String FEATURE_PANEL_NAME = "FeaturePanel";


    public TrackManager(IGV igv) {

        this.igv = igv;
        loader = new TrackLoader();
    }

    public Track getGeneTrack() {
        return geneTrack;
    }


    public SequenceTrack getSequenceTrack() {
        return sequenceTrack;
    }

    public Set<TrackType> getLoadedTypes() {
        return loadedTypes;
    }


    public void putScrollPane(String name, TrackPanelScrollPane sp) {
        trackPanelScrollPanes.put(name, sp);
    }

    public TrackPanelScrollPane getScrollPane(String name) {
        return trackPanelScrollPanes.get(name);
    }

    public Collection<TrackPanelScrollPane> getTrackPanelScrollPanes() {
        return trackPanelScrollPanes.values();
    }

    public void removeScrollPane(String name) {
        trackPanelScrollPanes.remove(name);
    }

    public void clearScrollPanes() {
        trackPanelScrollPanes.clear();
    }


    public void setGroupByAttribute(String attributeName) {
        groupByAttribute = attributeName;
        groupTracksByAttribute();
    }

    public void reset() {
        groupByAttribute = "";
        TrackPanelScrollPane tsp = trackPanelScrollPanes.get(DATA_PANEL_NAME);
        if (tsp != null) {
            tsp.getTrackPanel().reset();
        }
        loadedTypes.clear();
    }


    public void groupTracksByAttribute() {
        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel trackPanel = tsp.getTrackPanel();
            trackPanel.groupTracksByAttribute(groupByAttribute);
        }
        //TODO We should go through all current panels not just the data panel
        //TrackPanelScrollPane tsp = trackPanelScrollPanes.get(DATA_PANEL_NAME);
        //if (tsp != null) {
        //    tsp.getTrackPanel().groupTracksByAttribute(groupByAttribute);
        //}
    }


    public void chromosomeChanged(String chrName) {
        for (Track t : getAllTracks(true)) {
            t.chromosomeChanged(chrName);
        }
    }

    /**
     * A (hopefully) temporary solution to force SAM track reloads,  until we have a better
     * dependency scheme
     */
    public void reloadSAMTracks() {
        for (Track t : getAllTracks(false)) {
            if (t instanceof AlignmentTrack) {
                ((AlignmentTrack) t).clearCaches();
            }
        }
        IGV.getInstance().repaintDataPanels();
    }


    public void sortAlignmentTracks(AlignmentTrack.SortOption option) {
        for (Track t : getAllTracks(false)) {
            if (t instanceof AlignmentTrack) {
                for (ReferenceFrame frame : FrameManager.getFrames()) {
                    ((AlignmentTrack) t).sortRows(option, frame);
                }
            }
        }
    }

    public void sortAlignmentTracks(AlignmentTrack.SortOption option, double location) {
        for (Track t : getAllTracks(false)) {
            if (t instanceof AlignmentTrack) {
                for (ReferenceFrame frame : FrameManager.getFrames()) {
                    ((AlignmentTrack) t).sortRows(option, frame, location);
                }
            }
        }
    }

    public void packAlignmentTracks() {
        for (Track t : getAllTracks(false)) {
            if (t instanceof AlignmentTrack) {
                for (ReferenceFrame frame : FrameManager.getFrames()) {
                    ((AlignmentTrack) t).packAlignments(frame);
                }
            }
        }
    }

    public void collapseTracks() {
        for (Track t : getAllTracks(true)) {
            t.setDisplayMode(Track.DisplayMode.COLLAPSED);
        }
    }


    public void expandTracks() {
        for (Track t : getAllTracks(true)) {
            t.setDisplayMode(Track.DisplayMode.EXPANDED);
        }
    }

    public void collapseTrack(String trackName) {
        for (Track t : getAllTracks(true)) {
            if (t.getName().equals(trackName)) {
                t.setDisplayMode(Track.DisplayMode.COLLAPSED);
            }
        }
    }


    public void expandTrack(String trackName) {
        for (Track t : getAllTracks(true)) {
            if (t.getName().equals(trackName)) {
                t.setDisplayMode(Track.DisplayMode.EXPANDED);
            }
        }
    }


    public void loadResources(Collection<ResourceLocator> locators) {

        //Set<TrackPanel> changedPanels = new HashSet();

        log.info("Loading" + locators.size() + " resources.");
        MessageCollection messages = new MessageCollection();


        for (ResourceLocator locator : locators) {

            // If its a local file, check explicitly for existence (rather than rely on exception)
            if (locator.isLocal()) {
                File trackSetFile = new File(locator.getPath());
                if (!trackSetFile.exists()) {
                    messages.append("File not found: " + locator.getPath() + "\n");
                    continue;
                }
            }

            try {
                List<Track> tracks = load(locator);
                if (tracks.size() > 0) {
                    String path = locator.getPath();

                    // Get an appropriate panel.  If its a VCF file create a new panel if the number of genotypes
                    // is greater than 10
                    TrackPanel panel = getPanelFor(locator);
                    if (path.endsWith(".vcf") || path.endsWith(".vcf.gz") ||
                            path.endsWith(".vcf4") || path.endsWith(".vcf4.gz")) {
                        Track t = tracks.get(0);
                        if (t instanceof VCFTrack && ((VCFTrack) t).getAllSamples().size() > 10) {
                            String newPanelName = "Panel" + System.currentTimeMillis();
                            panel = igv.addDataPanel(newPanelName).getTrackPanel();
                        }
                    }
                    panel.addTracks(tracks);
                }
            } catch (Throwable e) {
                log.error("Error loading tracks", e);
                messages.append(e.getMessage());
            }
        }

        groupTracksByAttribute();
        resetOverlayTracks();

        //TODO Throw the message window higher up. There should be a throw when you are using the UI,
        // not when you are connecting through Socket

        if (!messages.isEmpty()) {
            for (String message : messages.getMessages()) {
                MessageUtils.showMessage(message);
            }
        }
    }

    /**
     * Load a resource (track or sample attribute file)
     */

    public List<Track> load(ResourceLocator locator) {

        try {
            List<Track> newTracks = loader.load(locator, igv);
            if (newTracks.size() > 0) {
                for (Track track : newTracks) {
                    String fn = locator.getPath();
                    int lastSlashIdx = fn.lastIndexOf("/");
                    if (lastSlashIdx < 0) {
                        lastSlashIdx = fn.lastIndexOf("\\");
                    }
                    if (lastSlashIdx > 0) {
                        fn = fn.substring(lastSlashIdx);
                    }
                    track.setAttributeValue("NAME", track.getName());
                    track.setAttributeValue("DATA FILE", fn);
                    track.setAttributeValue("DATA TYPE", track.getTrackType().toString());


                }
            }
            return newTracks;

        } catch (DataLoadException dle) {
            throw dle;
        } catch (Exception e) {
            log.error(e);
            throw new DataLoadException(e.getMessage(), locator.getPath());
        }

    }

    /**
     * @param type
     */
    public void addLoadedType(TrackType type) {
        loadedTypes.add(type);
    }

    /**
     * Load the data file into the specified panel.   Triggered via drag and drop.
     *
     * @param file
     * @param panel
     * @return
     */
    public void load(File file, TrackPanel panel) {
        ResourceLocator locator = new ResourceLocator(file.getAbsolutePath());
        List<Track> tracks = load(locator);
        panel.addTracks(tracks);
    }

    /**
     * Reset the overlay tracks collection.  Currently the only overlayable track
     * type is Mutation.  This method finds all mutation tracks and builds a map
     * of key -> mutatinon track,  where the key is the specified attribute value
     * for linking tracks for overlay.
     * <p/>
     * The method also resets all tracks overlay properties from the user
     * user preference value.  That should probably be done in another method.
     */
    public void resetOverlayTracks() {
        overlayTracksMap.clear();

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.OVERLAY_TRACKS_KEY)) {
            String overlayAttribute = IGV.getInstance().getSession().getOverlayAttribute();
            if (overlayAttribute != null) {
                for (Track track : getAllTracks(false)) {
                    if (track != null) {
                        if (track.getTrackType() == UIConstants.overlayTrackType) {
                            String value = track.getAttributeValue(overlayAttribute);

                            if (value != null) {
                                List<Track> trackList = overlayTracksMap.get(value);

                                if (trackList == null) {
                                    trackList = new ArrayList();
                                    overlayTracksMap.put(value, trackList);
                                }

                                trackList.add(track);
                            }
                        }
                    }
                }
            }
        }

        boolean displayOverlays = IGV.getInstance().getSession().getDisplayOverlayTracks();

        for (Track track : getAllTracks(false)) {
            if (track != null) {
                if (track.getTrackType() == UIConstants.overlayTrackType) {
                    track.setOverlayVisible(displayOverlays);
                }
            }

        }
    }


    /**
     * Method description
     *
     * @return
     */
    public int getVisibleTrackCount() {
        int count = 0;
        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel tsv = tsp.getTrackPanel();
            count += tsv.getVisibleTrackCount();

        }
        return count;
    }

    public List<Track> getAllTracks(boolean includeGeneTrack) {
        List<Track> allTracks = new ArrayList<Track>();

        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel tsv = tsp.getTrackPanel();
            allTracks.addAll(tsv.getTracks());
        }
        if ((geneTrack != null) && !includeGeneTrack) {
            allTracks.remove(geneTrack);
        }
        if ((sequenceTrack != null) && !includeGeneTrack) {
            allTracks.remove(sequenceTrack);
        }

        return allTracks;
    }

    public void clearSelections() {
        for (Track t : getAllTracks(true)) {
            if (t != null)
                t.setSelected(false);
        }

    }

    public void setTrackSelections(Set<Track> selectedTracks) {
        for (Track t : getAllTracks(true)) {
            if (selectedTracks.contains(t)) {
                t.setSelected(true);
            }
        }
    }

    public void shiftSelectTracks(Track track) {
        List<Track> allTracks = getAllTracks(true);
        int clickedTrackIndex = allTracks.indexOf(track);
        // Find another track that is already selected.  The semantics of this
        // are not well defined, so any track will do
        int otherIndex = clickedTrackIndex;
        for (int i = 0; i < allTracks.size(); i++) {
            if (allTracks.get(i).isSelected() && i != clickedTrackIndex) {
                otherIndex = i;
                break;
            }
        }

        int left = Math.min(otherIndex, clickedTrackIndex);
        int right = Math.max(otherIndex, clickedTrackIndex);
        for (int i = left; i <= right; i++) {
            allTracks.get(i).setSelected(true);
        }
    }

    public void toggleTrackSelections(Set<Track> selectedTracks) {
        for (Track t : getAllTracks(true)) {
            if (selectedTracks.contains(t)) {
                t.setSelected(!t.isSelected());
            }
        }
    }

    public Collection<Track> getSelectedTracks() {
        HashSet<Track> selectedTracks = new HashSet();
        for (Track t : getAllTracks(true)) {
            if (t != null && t.isSelected()) {
                selectedTracks.add(t);
            }
        }
        return selectedTracks;

    }

    /**
     * Return the complete set of unique DataResourcLocators currently loaded
     *
     * @return
     */
    public Set<ResourceLocator> getDataResourceLocators() {
        HashSet<ResourceLocator> locators = new HashSet();

        for (Track track : getAllTracks(false)) {
            ResourceLocator locator = track.getResourceLocator();

            if (locator != null) {
                locators.add(locator);
            }
        }

        return locators;

    }

    /**
     * Method description
     *
     * @param newHeight
     */
    public void setAllTrackHeights(int newHeight) {
        for (Track track : this.getAllTracks(false)) {
            track.setHeight(newHeight);
        }

    }


    public void unloadTracks(Collection<ResourceLocator> locators) {

        Collection<Track> tracksToRemove = new ArrayList();
        HashSet<ResourceLocator> locatorSet = new HashSet(locators);
        for (Track t : getAllTracks(true)) {
            if (locatorSet.contains(t.getResourceLocator())) {
                tracksToRemove.add(t);
            }
        }
        removeTracks(tracksToRemove);
    }

    /**
     * Method description
     *
     * @param tracksToRemove
     */
    public void removeTracks(Collection<Track> tracksToRemove) {

        // Make copy of list as we will be modifying the original in the loop
        ArrayList<TrackPanelScrollPane> panes = new ArrayList(getTrackPanelScrollPanes());
        for (TrackPanelScrollPane tsp : panes) {
            TrackPanel trackPanel = tsp.getTrackPanel();
            trackPanel.removeTracks(tracksToRemove);

            if (!trackPanel.hasTracks()) {
                igv.removeDataPanel(tsp.getTrackPanelName());
            }
        }

        for (Track t : tracksToRemove) {
            if (t instanceof DragListener) {
                DragEventManager.getInstance().removeDragListener((DragListener) t);
            }

        }
    }

    /**
     * Refresh the data in all tracks that respond to this message (note: most do not,  it was added early
     * on for copy number data and is kept for backward compatibility).
     */
    public void refreshData() {

        long t0 = System.currentTimeMillis();
        for (Track track : getAllTracks(false)) {
            track.refreshData(t0);
        }
    }

    /**
     * Sort all groups (data and feature) by attribute value(s).  Tracks are
     * sorted within groups.
     *
     * @param attributeNames
     * @param ascending
     */
    public void sortAllTracksByAttributes(final String attributeNames[], final boolean[] ascending) {
        assert attributeNames.length == ascending.length;

        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel tsv = tsp.getTrackPanel();
            tsv.sortTracksByAttributes(attributeNames, ascending);
        }
    }

    public void sortAllTracksByPosition(List<String> trackIds) {
        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel tsv = tsp.getTrackPanel();
            tsv.sortTracksByPosition(trackIds);
        }
    }

    /**
     * Return a DataPanel appropriate for the resource type
     *
     * @param locator
     * @return
     */
    public TrackPanel getPanelFor(ResourceLocator locator) {
        String path = locator.getPath().toLowerCase();
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            return igv.getDataPanel(DATA_PANEL_NAME);
        } else if (path.endsWith(".sam") || path.endsWith(".bam") ||
                path.endsWith(".sam.list") || path.endsWith(".bam.list") ||
                path.endsWith(".aligned") || path.endsWith(".sorted.txt")) {

            String newPanelName = "Panel" + System.currentTimeMillis();
            return igv.addDataPanel(newPanelName).getTrackPanel();
            //} else if (path.endsWith(".vcf") || path.endsWith(".vcf.gz") ||
            //        path.endsWith(".vcf4") || path.endsWith(".vcf4.gz")) {
            //    String newPanelName = "Panel" + System.currentTimeMillis();
            //    return igv.addDataPanel(newPanelName).getTrackPanel();
        } else {
            return getDefaultPanel(locator);
        }
    }


    private TrackPanel getDefaultPanel(ResourceLocator locator) {

        if (locator.getType() != null && locator.getType().equalsIgnoreCase("das")) {
            return igv.getDataPanel(FEATURE_PANEL_NAME);
        }

        String filename = locator.getPath().toLowerCase();

        if (filename.endsWith(".txt") || filename.endsWith(".tab") || filename.endsWith(
                ".xls") || filename.endsWith(".gz")) {
            filename = filename.substring(0, filename.lastIndexOf("."));
        }


        if (filename.contains("refflat") || filename.contains("ucscgene") ||
                filename.contains("genepred") || filename.contains("ensgene") ||
                filename.contains("refgene") ||
                filename.endsWith("gff") || filename.endsWith("gtf") ||
                filename.endsWith("gff3") || filename.endsWith("embl") ||
                filename.endsWith("bed") || filename.endsWith("gistic") ||
                filename.endsWith("bedz") || filename.endsWith("repmask") ||
                filename.contains("dranger")) {
            return igv.getDataPanel(FEATURE_PANEL_NAME);
        } else {
            return igv.getDataPanel(DATA_PANEL_NAME);
        }
    }

    /**
     * Sort all groups (data and feature) by a computed score over a region.  The
     * sort is done twice (1) groups are sorted with the featureGroup, and (2) the
     * groups themselves are sorted.
     *
     * @param region
     * @param type
     */
    public void sortByRegionScore(RegionOfInterest region,
                                  final RegionScoreType type,
                                  final ReferenceFrame frame) {

        final RegionOfInterest r = region == null ? new RegionOfInterest(frame.getChrName(), (int) frame.getOrigin(),
                (int) frame.getEnd() + 1, frame.getName()) : region;
        NamedRunnable runnable = new NamedRunnable() {

            public String getName() {
                return "Sort";
            }

            public void run() {
                for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
                    TrackPanel tsv = tsp.getTrackPanel();
                    tsv.sortByRegionsScore(r, type, frame);
                }
                IGV.getMainFrame().repaint();
            }
        };

        LongRunningTask.submit(runnable);

    }

    /**
     * Method description
     *
     * @return
     */
    public String getGroupByAttribute() {
        return groupByAttribute;
    }


    /**
     * Method description
     *
     * @param track
     * @return
     */
    public List<Track> getOverlayTracks(Track track) {
        String overlayAttribute = IGV.getInstance().getSession().getOverlayAttribute();
        String value = track.getAttributeValue(overlayAttribute);
        return overlayTracksMap.get(value);

    }

    /**
     * Method description
     *
     * @param genome
     */
    public void createGeneTrack(Genome genome) {

        FeatureTrack geneFeatureTrack = null;
        String genomeId = genome.getId();
        GenomeDescriptor genomeDescriptor = IGV.getInstance().getGenomeManager().getGenomeDescriptor(genomeId);
        AsciiLineReader reader = getGeneReader(genomeDescriptor);
        if (reader != null) {
            try {

                String geneFilename = genomeDescriptor.getGeneFileName();
                FeatureParser parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(geneFilename), genome);
                if (parser == null) {
                    MessageUtils.showMessage("ERROR: Unrecognized annotation file format: " + geneFilename +
                            "<br>Annotations for genome: " + genomeId + " will not be loaded.");
                } else {
                    List<org.broad.tribble.Feature> genes = parser.loadFeatures(reader);
                    String name = genomeDescriptor.getGeneTrackName();
                    if (name == null) name = "Genes";
                    geneFeatureTrack = new FeatureTrack(name, new FeatureCollectionSource(genes, genome));

                    geneFeatureTrack.setMinimumHeight(5);
                    geneFeatureTrack.setHeight(35);
                    geneFeatureTrack.setPreferredHeight(35);
                    geneFeatureTrack.setRendererClass(IGVFeatureRenderer.class);
                    geneFeatureTrack.setColor(Color.BLUE.darker());
                    TrackProperties props = parser.getTrackProperties();
                    if (props != null) {
                        geneFeatureTrack.setProperties(parser.getTrackProperties());
                    }
                    geneFeatureTrack.setUrl(genome.getAnnotationURL());
                }

            } catch (Exception e) {
                log.error("Error loading geneManager", e);
            }
            finally {

                if (reader != null) {
                    reader.close();
                }

            }
        }


        SequenceTrack seqTrack = new SequenceTrack("Reference sequence");
        if (geneFeatureTrack != null) {
            setGenomeTracks(geneFeatureTrack, seqTrack);
        } else {
            setGenomeTracks(null, seqTrack);
        }

    }

    /**
     * Replace current gene track with new one.  This is called upon switching genomes
     *
     * @param newGeneTrack
     * @param newSeqTrack
     */
    private void setGenomeTracks(Track newGeneTrack, SequenceTrack newSeqTrack) {

        boolean foundSeqTrack = false;
        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel tsv = tsp.getTrackPanel();
            foundSeqTrack = tsv.replaceTrack(sequenceTrack, newSeqTrack);
            if (foundSeqTrack) {
                break;
            }
        }

        boolean foundGeneTrack = false;
        for (TrackPanelScrollPane tsp : getTrackPanelScrollPanes()) {
            TrackPanel tsv = tsp.getTrackPanel();
            foundGeneTrack = tsv.replaceTrack(geneTrack, newGeneTrack);
            if (foundGeneTrack) {
                break;
            }
        }


        if (!foundGeneTrack || !foundSeqTrack) {
            TrackPanel panel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY) ?
                    igv.getDataPanel(DATA_PANEL_NAME) : igv.getDataPanel(FEATURE_PANEL_NAME);

            if (!foundSeqTrack) panel.addTrack(newSeqTrack);
            if (!foundGeneTrack && newGeneTrack != null) panel.addTrack(newGeneTrack);

        }

        // Keep a reference to this track so it can be removed
        geneTrack = newGeneTrack;
        sequenceTrack = newSeqTrack;

    }


    public static AsciiLineReader getGeneReader(GenomeDescriptor genomeDescriptor) {

        InputStream is = null;
        try {

            InputStream inputStream = genomeDescriptor.getGeneStream();
            if (inputStream == null) {
                return null;
            }

            AsciiLineReader reader = null;
            if (genomeDescriptor.isGeneFileGZipFormat()) {
                is = new GZIPInputStream(inputStream);
                reader = new AsciiLineReader(is);
            } else {
                is = new BufferedInputStream(inputStream);
                reader = new AsciiLineReader(is);
            }

            return reader;
        } catch (IOException ex) {
            log.warn("Error loading the genome!", ex);
            return null;
        }

    }


    /*
  public void updateTrackPositions(String panelName) {
  int regionY = 0;
  Collection<TrackGroup> groups = trackGroups.get(panelName).values();
  for (TrackGroup group : groups) {
  if (group.isVisible()) {
  if (groups.size() > 1) {
  regionY += UIStringConstants.groupGap;
  }
  int previousRegion = -1;
  for (Track track : group.getTracks()) {
  int trackHeight = track.getHeight();
  if (track.isVisible()) {
  if (isDragging) {
  Color currentColor = g.getColor();
  if (dragY >= previousRegion && dragY <= regionY) {
  g.setColor(Color.GRAY);
  g.drawLine(0, regionY, getWidth(), regionY);
  regionY++;
  g.setColor(currentColor);
  dropTarget = track;
  } else if (dragY >= (regionY + trackHeight)) {
  g.setColor(Color.GRAY);
  g.drawLine(0, regionY + trackHeight,
  getWidth(), regionY + trackHeight);
  trackHeight--;
  g.setColor(currentColor);
  dropTarget = null;
  }
  }
  previousRegion = regionY;
  if (regionY + trackHeight >= visibleRect.y) {
  int width = getWidth();
  int height = track.getHeight();
  Rectangle region = new Rectangle(regionX,
  regionY, width, height);
  addMousableRegion(new MouseableRegion(region, track));
  draw(graphics2D, track, regionX, regionY, width,
  height, visibleRect);
  }
  regionY += trackHeight;
  }
  }
  if (group.isDrawBorder()) {
  g.drawLine(0, regionY, getWidth(), regionY);
  }
  }
  }
  }
   * */
    /**
     * Method description
     *
     */
    /*
    public void dumpData() {
    PrintWriter pw = null;
    try {
    String chr = ReferenceFrame.getInstance().getChromosome().getName();
    double x = ReferenceFrame.getInstance().getOrigin();
    pw = new PrintWriter(new FileWriter("DataDump.tab"));
    pw.println("Sample\tCopy Number\tExpression\tMethylation");
    for (String name : groupsMap.keySet()) {
    pw.print(name);
    float cn = 0;
    float exp = 0;
    float meth = 0;
    int nCn = 0;
    int nExp = 0;
    int nMeth = 0;
    for (Track t : groupsMap.get(name).getTracks()) {
    LocusScore score = ((DataTrack) t).getLocusScoreAt(chr, x);
    if ((score != null) && !Float.isNaN(score.getScore())) {
    if (t.getTrackType() == TrackType.COPY_NUMBER) {
    nCn++;
    cn += score.getScore();
    } else if (t.getTrackType() == TrackType.GENE_EXPRESSION) {
    nExp++;
    exp += score.getScore();
    } else if (t.getTrackType() == TrackType.DNA_METHYLATION) {
    nMeth++;
    meth += score.getScore();
    }
    }
    }
    pw.print("\t");
    if (nCn > 0) {
    pw.print(cn / nCn);
    }
    pw.print("\t");
    if (nExp > 0) {
    pw.print(exp / nExp);
    }
    pw.print("\t");
    if (nMeth > 0) {
    pw.print(cn / meth);
    }
    pw.println();
    }
    for (Track track : getTracksForPanel(false)) {
    if (track instanceof DataTrack) {
    LocusScore score = ((DataTrack) track).getLocusScoreAt(chr,
    x);
    if (score != null) {
    pw.println(
    track.getName() + "\t" + track.getTrackType() + "\t" + score.getScore());
    }
    }
    }
    } catch (IOException ex) {
    ex.printStackTrace();
    } finally {
    pw.close();
    }
    }     * */

}
