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
 * GenomeManager.java
 *
 * Created on November 9, 2007, 9:12 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;


import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.load.*;
import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ucsc.hub.Hub;
import org.broad.igv.ucsc.hub.HubParser;
import org.broad.igv.ucsc.hub.TrackSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.PanelName;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.genome.GenomeListManager;
import org.broad.igv.ui.genome.GenomeListItem;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.net.SocketException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author jrobinso
 */
public class GenomeManager {

    public static final String SELECT_ANNOTATIONS_MESSAGE = "Select default annotation tracks for this genome.  " +
            "You can change these selections later using the 'Genomes > Select Genome Annotations...' menu.";
    private static Logger log = LogManager.getLogger(GenomeManager.class);

    private static GenomeManager theInstance;

    private static GenomeListManager genomeListManager;

    private Genome currentGenome;


    /**
     * Map from genomeID -> GenomeTableRecord
     * ID comparison will be case insensitive
     */

    public synchronized static GenomeManager getInstance() {
        if (theInstance == null) {
            theInstance = new GenomeManager();
        }
        return theInstance;
    }

    private GenomeManager() {
        genomeListManager = GenomeListManager.getInstance();
        GenomeLoader.loadSequenceMap();
    }


    public String getGenomeId() {
        return currentGenome == null ? null : currentGenome.getId();
    }

    /**
     * IGV always has exactly 1 genome loaded at a time.
     * This returns the currently loaded genome
     *
     * @return
     * @api
     */
    public Genome getCurrentGenome() {
        return currentGenome;
    }

    /**
     * Load a genome by ID, which might be a file path or URL
     *
     * @param genomeId - ID for an IGV hosted genome, or file path or url
     * @return boolean flag indicating success
     * @throws IOException
     */
    public boolean loadGenomeById(String genomeId) throws IOException {

        final Genome currentGenome = getCurrentGenome();
        if (currentGenome != null && genomeId.equals(currentGenome.getId())) {
            return false;
        }

        String genomePath;
        if (org.broad.igv.util.ParsingUtils.fileExists(genomeId)) {
            genomePath = genomeId;
        } else {
            GenomeListItem item = getGenomeTableRecord(genomeId);
            if (item == null) {
                MessageUtils.showMessage("Could not locate genome with ID: " + genomeId);
                return false;
            } else {
                genomePath = item.getPath();
            }
        }
        return loadGenome(genomePath) != null; // monitor[0]);
    }

    /**
     * The main load method -- loads a genome from a file or url path.  Note this is a long running operation and
     * should not be done on the Swing event thread as it will block the UI.
     * <p>
     * NOTE: The member 'currentGenome' is set here as a side effect.
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    public Genome loadGenome(String genomePath) throws IOException {

        WaitCursorManager.CursorToken cursorToken = null;
        try {
            log.info("Loading genome: " + genomePath);
            if (IGV.hasInstance()) {
                IGVMenuBar.getInstance().disableTracksMenu();
                IGV.getInstance().setStatusBarMessage("<html><font color=blue>Loading genome</font></html>");
                cursorToken = WaitCursorManager.showWaitCursor();
            }

            Genome newGenome = GenomeLoader.getLoader(genomePath).loadGenome();

            // Load user-defined chr aliases, if any.  This is done last so they have priority
            final File genomeCacheDirectory = DirectoryManager.getGenomeCacheDirectory();
            try {
                String aliasPath = (new File(genomeCacheDirectory, newGenome.getId() + "_alias.tab")).getAbsolutePath();
                if (!(new File(aliasPath)).exists()) {
                    aliasPath = (new File(genomeCacheDirectory, newGenome.getId() + "_alias.tab.txt")).getAbsolutePath();
                }
                if ((new File(aliasPath)).exists()) {
                    newGenome.addChrAliases(ChromAliasParser.loadChrAliases(aliasPath));
                }
            } catch (Exception e) {
                log.error("Failed to load user defined alias", e);
            }


            if (IGV.hasInstance()) {
                IGV.getInstance().resetSession(null);
            }

            // Add an entry to the pulldown
            GenomeListItem genomeListItem = new GenomeListItem(newGenome.getDisplayName(), genomePath, newGenome.getId());
            GenomeListManager.getInstance().addGenomeItem(genomeListItem);

            setCurrentGenome(newGenome);

            return currentGenome;

        } catch (SocketException e) {
            throw new RuntimeException("Server connection error", e);
        } finally {
            if (IGV.hasInstance()) {
                IGV.getInstance().setStatusBarMessage("");
                WaitCursorManager.removeWaitCursor(cursorToken);
                IGVMenuBar.getInstance().enableTracksMenu();
            }
        }
    }

    public void setCurrentGenome(Genome newGenome) {

        this.currentGenome = newGenome;

        // hasInstance() check to filters unit test
        if (IGV.hasInstance()) {
            IGV.getInstance().goToLocus(newGenome.getHomeChromosome()); //  newGenome.getDefaultPos());
            FrameManager.getDefaultFrame().setChromosomeName(newGenome.getHomeChromosome(), true);

            restoreGenomeTracks(newGenome);

            IGV.getInstance().resetFrames();
            IGV.getInstance().getSession().clearHistory();

            if (newGenome != Genome.nullGenome()) {
                // This should only occur on startup failure
                PreferencesManager.getPreferences().setLastGenome(newGenome.getId());
            }

            if (PreferencesManager.getPreferences().getAsBoolean(Constants.CIRC_VIEW_ENABLED) && CircularViewUtilities.ping()) {
                CircularViewUtilities.changeGenome(newGenome);
            }

            IGVEventBus.getInstance().post(new GenomeChangeEvent(newGenome));
        }
    }

    /**
     * @param genome
     */
    public void restoreGenomeTracks(Genome genome) {

        IGV.getInstance().setSequenceTrack();

        // Fetch the gene track, defined by .genome files.  In this format the genome data is encoded in the .genome file
        FeatureTrack geneFeatureTrack = genome.getGeneTrack();   // Used for .genome and .gbk formats.  Otherwise null
        if (geneFeatureTrack != null) {
            IGV.getInstance().addTrack(geneFeatureTrack, PanelName.ANNOTATION_PANEL.getName());
        }

        List<ResourceLocator> resources = genome.getAnnotationResources();
        List<Track> annotationTracks = new ArrayList<>();
        if (resources != null) {
            for (ResourceLocator locator : resources) {
                try {
                    if(locator != null) {
                        locator.setPanelName(PanelName.ANNOTATION_PANEL.getName());
                        List<Track> tracks = IGV.getInstance().load(locator);
                        annotationTracks.addAll(tracks);
                    }
                } catch (DataLoadException e) {
                    log.error("Error loading genome annotations", e);
                }
            }
        }

        if (annotationTracks.size() > 0) {
            IGV.getInstance().addTracks(annotationTracks);
            for (Track track : annotationTracks) {
                ResourceLocator locator = track.getResourceLocator();
                if(locator != null) {
                    String fn = "";
                    if (locator != null) {
                        fn = locator.getPath();
                        int lastSlashIdx = fn.lastIndexOf("/");
                        if (lastSlashIdx < 0) {
                            lastSlashIdx = fn.lastIndexOf("\\");
                        }
                        if (lastSlashIdx > 0) {
                            fn = fn.substring(lastSlashIdx + 1);
                        }
                    }
                }
            }
        }

        IGV.getInstance().revalidateTrackPanels();
    }

    /**
     * Update the annotation tracks for the current genome.  This will prompt the user to select tracks from the
     * genomes default hub.  The selected tracks will be saved in the genome config file, and loaded.  Deselected
     * tracks will be removed.
     *
     * @throws IOException
     */
    public void updateAnnotations() throws IOException {
        if (currentGenome != null) {
            GenomeConfig config = currentGenome.getConfig();

            if (config != null) {

                final List<TrackConfig> trackConfigs = config.getTrackConfigs();
                List<String> currentAnnotationPaths = trackConfigs == null ? Collections.EMPTY_LIST :
                        trackConfigs.stream().map(t -> t.url).toList();

                String message = "Select defaul tannoations for " + config.getName();
                List<TrackConfig> selectedConfigs = selectAnnotationTracks(config, message);
                if (selectedConfigs == null) {
                    return;
                }

                config.setTracks(selectedConfigs);
                GenomeDownloadUtils.saveLocalGenome(config);

                Set<String> selectedTrackPaths = selectedConfigs.stream().map(t -> t.url).collect(Collectors.toSet());

                // Unload deselected tracks
                Set<String> pathsToRemove = new HashSet<>();
                for (String p : currentAnnotationPaths) {
                    if (!selectedTrackPaths.contains(p)) {
                        pathsToRemove.add(p);
                    }
                }
                IGV.getInstance().deleteTracksByPath(pathsToRemove);

                // Load selected tracks.Filter out tracks already loaded
                Set<String> loadedTrackPaths = IGV.getInstance().getAllTracks().stream()
                        .filter(t -> t.getResourceLocator() != null)
                        .map(t -> t.getResourceLocator().getPath())
                        .collect(Collectors.toSet());
                List<TrackConfig> tracksToLoad = selectedConfigs.stream()
                        .filter(trackConfig -> !loadedTrackPaths.contains(trackConfig.url))
                        .collect(Collectors.toList());

                List<ResourceLocator> locators = tracksToLoad.stream().map(t -> ResourceLocator.fromTrackConfig(t)).toList();
                for (ResourceLocator locator : locators) {
                    locator.setPanelName(PanelName.ANNOTATION_PANEL.getName());
                }

                IGV.getInstance().loadTracks(locators);
            }
        }
    }

    /**
     * Prompt the user to select annotation tracks from the genome's default hub.
     *
     * @param config
     * @return
     * @throws IOException
     */

    public static List<TrackConfig> selectAnnotationTracks(GenomeConfig config, String message) throws IOException {

        String annotationHub = config.getHubs().get(0);  // IGV convention
        Hub hub = HubParser.loadHub(annotationHub);

        Set<String> currentSelections = config.getTrackConfigs() == null ? Collections.emptySet() :
                config.getTrackConfigs().stream()
                        .map(trackConfig -> trackConfig.url)
                        .collect(Collectors.toSet());
        TrackSelectionDialog dlg = TrackSelectionDialog.getTrackHubSelectionDialog(hub, config.getUcscID(), currentSelections, true, message);
        try {
            UIUtilities.invokeAndWaitOnEventThread(() -> dlg.setVisible(true));
            if (dlg.isCanceled()) {
                return null;
            } else {
                return dlg.getSelectedConfigs();
            }
        } catch (Exception e) {
            log.error("Error opening or using TrackHubSelectionDialog: " + e.getMessage());
            return null;
        }
    }

    /**
     * Delete the specified genome files
     *
     * @param removedValuesList
     */
    public void deleteDownloadedGenomes(List<GenomeListItem> removedValuesList) throws IOException {

        for (GenomeListItem item : removedValuesList) {

            String loc = item.getPath();
            File genomeFile = new File(loc);

            if (genomeFile.exists() && DirectoryManager.isChildOf(DirectoryManager.getGenomeCacheDirectory(), genomeFile)) {

                genomeFile.delete();

                // Delete associated data files
                File dataFileDirectory = new File(DirectoryManager.getGenomeCacheDirectory(), item.getId());
                File localFasta = DotGenomeUtils.getLocalFasta(item.getId());  //  (Legacy .genome convention)

                if ((dataFileDirectory.isDirectory() || localFasta != null) &&
                        MessageUtils.confirm("Delete downloaded data files?")) {

                    if (dataFileDirectory.isDirectory()) {
                        try (Stream<Path> paths = Files.walk(dataFileDirectory.toPath())) {
                            paths.sorted(Comparator.reverseOrder())
                                    .map(Path::toFile)
                                    .forEach(File::delete);
                        }
                        dataFileDirectory.delete();
                    }

                    if (localFasta != null) {
                        // If fasta file is in the "igv/genomes" directory delete it
                        DotGenomeUtils.removeLocalFasta(item.getId());
                        if (DirectoryManager.isChildOf(DirectoryManager.getGenomeCacheDirectory(), localFasta)) {
                            if (MessageUtils.confirm("Delete fasta file: " + localFasta.getAbsolutePath() + "?")) {
                                localFasta.delete();
                                File indexFile = new File(localFasta.getAbsolutePath() + ".fai");
                                if (indexFile.exists()) {
                                    indexFile.delete();
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    /**
     * Searches through currently loaded GenomeTableRecords and returns
     * that with a matching ID. If not found, searches server and
     * user defined lists
     *
     * @param genomeId
     * @return
     */
    public GenomeListItem getGenomeTableRecord(String genomeId) {

        GenomeListItem matchingItem = GenomeListManager.getInstance().getGenomeTableRecord(genomeId);
        if (matchingItem == null) {
            // If genome archive was not found, search hosted genomes
            matchingItem = HostedGenomes.getGenomeListItem(genomeId);
        }
        return matchingItem;
    }

    // Setter provided for unit tests
    public void setCurrentGenomeForTest(Genome genome) {
        this.currentGenome = genome;
    }


}
