/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

/*
 * IGV.java
 *
 * Represents an IGV instance.
 *
 * Note:  Currently, only one instance is allowed per JVM.
 *
 */
package org.broad.igv.ui;

import com.jidesoft.swing.JideSplitPane;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.*;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.lists.Preloader;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.peaks.PeakCommandBar;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionReader;
import org.broad.igv.track.*;

import static org.broad.igv.ui.WaitCursorManager.CursorToken;

import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.event.AlignmentTrackEventListener;
import org.broad.igv.ui.event.TrackGroupEvent;
import org.broad.igv.ui.event.TrackGroupEventListener;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;

import static org.broad.igv.ui.util.SnapshotUtilities.*;

import org.broad.igv.ui.util.ProgressMonitor;


import org.broad.igv.util.*;
import org.broad.igv.variant.VariantTrack;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.SeekableFileStream;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.*;
import java.lang.ref.SoftReference;
import java.net.*;
import java.util.*;
import java.util.List;

/**
 * Represents an IGV instance, consisting of a main window and associated model.
 *
 * @author jrobinso
 */
public class IGV {

    private static Logger log = Logger.getLogger(IGV.class);
    private static IGV theInstance;
    static List<IGV> instances = new LinkedList();

    // Window components
    private Frame mainFrame;
    private JRootPane rootPane;
    private IGVContentPane contentPane;
    private IGVMenuBar menuBar;

    private StatusWindow statusWindow;

    // Glass panes
    Component glassPane;
    GhostGlassPane dNdGlassPane;

    // Cursors
    public static Cursor fistCursor;
    public static Cursor zoomInCursor;
    public static Cursor zoomOutCursor;
    public static Cursor dragNDropCursor;

    //Session session;
    Session session;

    private GenomeManager genomeManager;
    /**
     * The gene track for the current genome, rendered in the FeaturePanel
     */
    private Track geneTrack;

    /**
     * The sequence track for the current genome
     */
    private SequenceTrack sequenceTrack;

    /**
     * Attribute used to group tracks.  Normally "null".  Set from the "Tracks" menu.
     */
    private String groupByAttribute = null;


    private Map<String, List<Track>> overlayTracksMap = new HashMap();
    private Set<Track> overlaidTracks = new HashSet();

    public static final String DATA_PANEL_NAME = "DataPanel";
    public static final String FEATURE_PANEL_NAME = "FeaturePanel";


    // Misc state
    private LinkedList<String> recentSessionList = new LinkedList<String>();
    private boolean isExportingSnapshot = false;
    private boolean startupComplete = false;


    /**
     *
     */
    Collection<SoftReference<TrackGroupEventListener>> groupListeners =
            Collections.synchronizedCollection(new ArrayList<SoftReference<TrackGroupEventListener>>());

    Collection<SoftReference<AlignmentTrackEventListener>> alignmentTrackListeners =
            Collections.synchronizedCollection(new ArrayList<SoftReference<AlignmentTrackEventListener>>());


    public static IGV createInstance(Frame frame) {
        if (theInstance != null) {
            throw new RuntimeException("Only a single instance is allowed.");
        }
        theInstance = new IGV(frame);
        return theInstance;
    }

    public static IGV getInstance() {
        if (theInstance == null) {
            throw new RuntimeException("IGV has not been initialized.  Must call createInstance(Frame) first");
        }
        return theInstance;
    }


    public static IGV getFirstInstance() {
        if (instances.isEmpty()) {
            throw new RuntimeException("IGV has not been initialized.  Must call createInstance(Frame) first");
        }
        return instances.get(0);
    }


    public static boolean hasInstance() {
        return theInstance != null;
    }


    public static JRootPane getRootPane() {
        return getInstance().rootPane;
    }

    public static Frame getMainFrame() {
        return getInstance().mainFrame;
    }


    /**
     * Creates new form IGV
     */
    private IGV(Frame frame) {

        theInstance = this;
        instances.add(this);

        genomeManager = new GenomeManager(this);

        mainFrame = frame;
        mainFrame.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent windowEvent) {
                windowCloseEvent();
            }

            @Override
            public void windowClosed(WindowEvent windowEvent) {
                windowCloseEvent();
            }

            private void windowCloseEvent() {
                instances.remove(this);
                PreferenceManager.getInstance().setApplicationFrameBounds(rootPane.getBounds());

            }

            @Override
            public void windowLostFocus(WindowEvent windowEvent) {
                // Start & stop tooltip manager to force any tooltip windows to close.
                ToolTipManager.sharedInstance().setEnabled(false);
                ToolTipManager.sharedInstance().setEnabled(true);
                IGVPopupMenu.closeAll();
            }


            @Override
            public void windowDeactivated(WindowEvent windowEvent) {
                // Start & stop tooltip manager to force any tooltip windows to close.
                ToolTipManager.sharedInstance().setEnabled(false);
                ToolTipManager.sharedInstance().setEnabled(true);
                IGVPopupMenu.closeAll();
            }

            @Override
            public void windowActivated(WindowEvent windowEvent) {

            }

            @Override
            public void windowGainedFocus(WindowEvent windowEvent) {

            }
        });


        session = new Session(null);

        // Create cursors
        createHandCursor();
        createZoomCursors();
        createDragAndDropCursor();

        // Create components
        mainFrame.setTitle(UIConstants.APPLICATION_NAME);

        if (mainFrame instanceof JFrame) {
            JFrame jf = (JFrame) mainFrame;
            rootPane = jf.getRootPane();
        } else {
            rootPane = new JRootPane();
            mainFrame.add(rootPane);

        }
        contentPane = new IGVContentPane(this);
        menuBar = new IGVMenuBar();

        rootPane.setContentPane(contentPane);
        rootPane.setJMenuBar(menuBar);
        glassPane = rootPane.getGlassPane();
        glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        glassPane.addMouseListener(new MouseAdapter() {
        });
        dNdGlassPane = new GhostGlassPane();

        mainFrame.pack();

        // Set the application's previous location and size
        Dimension screenBounds = Toolkit.getDefaultToolkit().getScreenSize();
        Rectangle applicationBounds = PreferenceManager.getInstance().getApplicationFrameBounds();
        if (applicationBounds == null || applicationBounds.getMaxX() > screenBounds.getWidth() ||
                applicationBounds.getMaxY() > screenBounds.getHeight()) {
            int width = Math.min(1150, (int) screenBounds.getWidth());
            int height = Math.min(800, (int) screenBounds.getHeight());
            applicationBounds = new Rectangle(0, 0, width, height);
        }
        mainFrame.setBounds(applicationBounds);

    }


    public void repaint() {
        mainFrame.repaint();
    }


    public GhostGlassPane getDnDGlassPane() {
        return dNdGlassPane;
    }

    public void startDnD() {
        rootPane.setGlassPane(dNdGlassPane);
        dNdGlassPane.setVisible(true);
    }

    public void endDnD() {
        rootPane.setGlassPane(glassPane);
        glassPane.setVisible(false);
    }

    public Dimension getPreferredSize() {
        return UIConstants.preferredSize;
    }


    public void addRegionOfInterest(RegionOfInterest roi) {
        session.addRegionOfInterestWithNoListeners(roi);
        RegionOfInterestPanel.setSelectedRegion(roi);
        doRefresh();
    }

    void beginROI(JButton button) {
        for (TrackPanel tp : getTrackPanels()) {
            TrackPanelScrollPane tsv = tp.getScrollPane();
            DataPanelContainer dpc = tsv.getDataPanel();
            for (Component c : dpc.getComponents()) {
                if (c instanceof DataPanel) {
                    DataPanel dp = (DataPanel) c;
                    RegionOfInterestTool regionOfInterestTool = new RegionOfInterestTool(dp, button);
                    dp.setCurrentTool(regionOfInterestTool);
                }
            }
        }


    }

    public void endROI() {
        for (TrackPanel tp : getTrackPanels()) {
            DataPanelContainer dp = tp.getScrollPane().getDataPanel();
            dp.setCurrentTool(null);
        }

    }


    public void chromosomeChangeEvent(String chrName) {
        chromosomeChangeEvent(chrName, true);
    }

    public void chromosomeChangeEvent(String chrName, boolean updateCommandBar) {

        contentPane.chromosomeChanged(chrName);
        repaintDataAndHeaderPanels(updateCommandBar);

    }

    /**
     * Repaint panels containing data, specifically the dataTrackPanel,
     * featureTrackPanel, and headerPanel.
     */
    public void repaintDataAndHeaderPanels() {
        repaintDataAndHeaderPanels(true);
    }

    /**
     * Repaint the header and data panels.
     * <p/>
     * Note:  If running in Batch mode a monitor is used to force synchrnous painting.  This is neccessary as the
     * paint() command triggers loading of data.  If allowed to proceed asynchronously the "snapshot" batch command
     * might execute before the data from a previous command has loaded.
     *
     * @param updateCommandBar
     */
    public void repaintDataAndHeaderPanels(boolean updateCommandBar) {
        if (Globals.isBatch()) {
            if (SwingUtilities.isEventDispatchThread()) {
                rootPane.paintImmediately(rootPane.getBounds());
            } else {
                synchronized (this) {
                    Runnable r = new Runnable() {
                        public void run() {
                            synchronized (IGV.this) {
                                rootPane.paintImmediately(rootPane.getBounds());
                                IGV.this.notify();
                            }
                        }
                    };
                    UIUtilities.invokeOnEventThread(r);
                    try {
                        // Wait a maximum of 5 minutes
                        this.wait(5 * 60 * 1000);
                    } catch (InterruptedException e) {
                        // Just continue
                    }
                }
            }
        } else {
            rootPane.repaint();
        }
        if (updateCommandBar) {
            contentPane.updateCurrentCoordinates();
        }
    }

    /**
     * Repaint the data panels.  Deprecated, but kept for backwards compatibility.
     */
    public void repaintDataPanels() {
        repaintDataAndHeaderPanels(false);
    }

    public void repaintNamePanels() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().getNamePanel().repaint();
        }

    }

    public void repaintStatusAndZoomSlider() {
        contentPane.getCommandBar().repaint();
    }


    public void selectGenomeFromList(String genome) {
        try {
            contentPane.getCommandBar().selectGenomeFromList(genome);
        } catch (FileNotFoundException e) {
            log.error("File not found while intializing genome!", e);
        } catch (NoRouteToHostException e) {
            log.error("Error while intializing genome!", e);
        }

    }


    public void doDefineGenome(ProgressMonitor monitor) {

        ProgressBar bar = null;
        File archiveFile = null;

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            GenomeBuilderDialog genomeBuilderDialog = new GenomeBuilderDialog(this, true);

            genomeBuilderDialog.setVisible(true);
            if (genomeBuilderDialog.isCanceled()) {
                return;
            }

            if (monitor != null) {
                bar = ProgressBar.showProgressDialog(mainFrame, "Defining Genome...", monitor, false);
            }

            String genomeZipLocation = genomeBuilderDialog.getGenomeArchiveLocation();
            String cytobandFileName = genomeBuilderDialog.getCytobandFileName();
            String refFlatFileName = genomeBuilderDialog.getRefFlatFileName();
            String fastaFileName = genomeBuilderDialog.getFastaFileName();
            String chrAliasFile = genomeBuilderDialog.getChrAliasFileName();
            String relativeSequenceLocation = genomeBuilderDialog.getSequenceLocation();
            String seqLocationOverride = genomeBuilderDialog.getSequenceLocationOverride();
            String genomeDisplayName = genomeBuilderDialog.getGenomeDisplayName();
            String genomeId = genomeBuilderDialog.getGenomeId();
            String genomeFileName = genomeBuilderDialog.getArchiveFileName();

            GenomeListItem genomeListItem = IGV.getInstance().getGenomeManager().defineGenome(
                    genomeZipLocation, cytobandFileName, refFlatFileName,
                    fastaFileName, chrAliasFile, relativeSequenceLocation, genomeDisplayName,
                    genomeId, genomeFileName, monitor, seqLocationOverride);

            if (genomeListItem != null) {
                enableRemoveGenomes();

                contentPane.getCommandBar().addToUserDefinedGenomeItemList(genomeListItem);
                contentPane.getCommandBar().selectGenomeFromListWithNoImport(genomeListItem.getId());
            }
            if (monitor != null) {
                monitor.fireProgressChange(100);
            }

        } catch (MaximumContigGenomeException e) {

            String genomePath = "";
            if (archiveFile != null) {
                genomePath = archiveFile.getAbsolutePath();
            }

            log.error("Failed to define genome: " + genomePath, e);

            JOptionPane.showMessageDialog(mainFrame, "Failed to define the current genome " +
                    genomePath + "\n" + e.getMessage());
        } catch (Exception e) {
            String genomePath = "";
            if (archiveFile != null) {
                genomePath = archiveFile.getAbsolutePath();
            }

            log.error("Failed to define genome: " + genomePath, e);
            MessageUtils.showMessage("Unexpected while importing a genome: " + e.getMessage());
        } finally {
            if (bar != null) {
                bar.close();
            }
            WaitCursorManager.removeWaitCursor(token);
        }
    }

    public GenomeListItem getGenomeSelectedInDropdown() {
        return contentPane.getCommandBar().getGenomeSelectedInDropdown();
    }

    /**
     * Gets the collection of genome display names currently in use.
     *
     * @return Set of display names.
     */
    public Collection<String> getGenomeDisplayNames() {
        return contentPane.getCommandBar().getGenomeDisplayNames();
    }

    public Collection<String> getGenomeIds() {
        return contentPane.getCommandBar().getGenomeIds();
    }


    /**
     * Load a .genome file directly.  This method really belongs in IGVMenuBar.
     *
     * @param monitor
     * @return
     */

    public void doLoadGenome(ProgressMonitor monitor) {

        ProgressBar bar = null;
        File file = null;
        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            File importDirectory = PreferenceManager.getInstance().getLastGenomeImportDirectory();
            if (importDirectory == null) {
                PreferenceManager.getInstance().setLastGenomeImportDirectory(Globals.getUserDirectory());
            }

            // Display the dialog
            file = FileDialogUtils.chooseFile("Load Genome", importDirectory, FileDialog.LOAD);

            // If a file selection was made
            if (file != null) {
                if (monitor != null) {
                    bar = ProgressBar.showProgressDialog(mainFrame, "Loading Genome...", monitor, false);
                }

                loadGenome(file.getAbsolutePath(), monitor);

            }
        } catch (Exception e) {
            MessageUtils.showMessage("<html>Error loading: " + file.getAbsolutePath() + "<br>" + e.getMessage());
            log.error("Error loading: " + file.getAbsolutePath(), e);
        } finally {
            WaitCursorManager.removeWaitCursor(token);
            if (monitor != null) {
                monitor.fireProgressChange(100);
            }

            if (bar != null) {
                bar.close();
            }
        }

    }

    public void loadGenome(String path, ProgressMonitor monitor) throws IOException {

        File file = new File(path);
        if (file.exists()) {
            File directory = file.getParentFile();
            PreferenceManager.getInstance().setLastGenomeImportDirectory(directory);
        }

        Genome genome = getGenomeManager().loadGenome(path, monitor);
        final String name = genome.getDisplayName();
        final String id = genome.getId();

        GenomeListItem genomeListItem = new GenomeListItem(name, path, id, true);
        getGenomeManager().addUserDefineGenomeItem(genomeListItem);

        contentPane.getCommandBar().addToUserDefinedGenomeItemList(genomeListItem);
        contentPane.getCommandBar().selectGenomeFromListWithNoImport(genomeListItem.getId());

    }


    public void enableExtrasMenu() {

        menuBar.enableExtrasMenu();
    }

    /**
     * Load a collection of tracks in a background thread.
     * <p/>
     * Note: Most of the code here is to adjust the scrollbars and split pane after loading
     *
     * @param locators
     */
    public void loadTracks(final Collection<ResourceLocator> locators) {

        contentPane.getStatusBar().setMessage("Loading ...");

        log.debug("Run loadTracks");


        if (locators != null && !locators.isEmpty()) {

            // NOTE:  this work CANNOT be done on the dispatch thread, it will potentially cause deadlock if
            // dialogs are opened or other Swing tasks are done.

            NamedRunnable runnable = new NamedRunnable() {
                public void run() {
                    // get current track count per panel.  Needed to detect which panels
                    // changed.  Also record panel sizes
                    final HashMap<TrackPanelScrollPane, Integer> trackCountMap = new HashMap();
                    final HashMap<TrackPanelScrollPane, Integer> panelSizeMap = new HashMap();
                    for (TrackPanel tp : getTrackPanels()) {
                        TrackPanelScrollPane sp = tp.getScrollPane();
                        trackCountMap.put(sp, sp.getDataPanel().getAllTracks().size());
                        panelSizeMap.put(sp, sp.getDataPanel().getHeight());
                    }

                    loadResources(locators);

                    double totalHeight = 0;
                    for (TrackPanel tp : getTrackPanels()) {
                        TrackPanelScrollPane sp = tp.getScrollPane();
                        if (trackCountMap.containsKey(sp)) {
                            int prevTrackCount = trackCountMap.get(sp).intValue();
                            if (prevTrackCount != sp.getDataPanel().getAllTracks().size()) {
                                int scrollPosition = panelSizeMap.get(sp);
                                if (prevTrackCount != 0 && sp.getVerticalScrollBar().isShowing()) {
                                    sp.getVerticalScrollBar().setMaximum(sp.getDataPanel().getHeight());
                                    sp.getVerticalScrollBar().setValue(scrollPosition);
                                }
                            }
                        }
                        // Give a maximum "weight" of 300 pixels to each panel.  If there are no tracks, give zero
                        if (sp.getTrackPanel().getTracks().size() > 0)
                            totalHeight += Math.min(300, sp.getTrackPanel().getPreferredPanelHeight());
                    }

                    // Adjust dividers for data panel.  The data panel divider can be
                    // zero if there are no data tracks loaded.
                    final JideSplitPane centerSplitPane = contentPane.getMainPanel().getCenterSplitPane();
                    int htotal = centerSplitPane.getHeight();
                    int y = 0;
                    int i = 0;
                    for (Component c : centerSplitPane.getComponents()) {
                        if (c instanceof TrackPanelScrollPane) {
                            final TrackPanel trackPanel = ((TrackPanelScrollPane) c).getTrackPanel();
                            if (trackPanel.getTracks().size() > 0) {
                                int panelWeight = Math.min(300, trackPanel.getPreferredPanelHeight());
                                int dh = (int) ((panelWeight / totalHeight) * htotal);
                                y += dh;
                            }
                            centerSplitPane.setDividerLocation(i, y);
                            i++;
                        }
                    }

                    SwingUtilities.invokeLater(new Runnable() {
                        public void run() {
                            contentPane.getMainPanel().doLayout();
                            contentPane.getStatusBar().setMessage("");
                        }
                    });
                    IGV.getInstance().showLoadedTrackCount();
                }

                public String getName() {
                    return "Load Tracks";
                }
            };

            LongRunningTask.submit(runnable);
        }

        log.debug("Finish loadTracks");

    }


    public void setGeneList(String listID) {
        setGeneList(listID, true);
    }

    public void setGeneList(final String listID, final boolean recordHistory) {

        //LongRunningTask.submit(new NamedRunnable() {
        //    public String getName() {
        //        return "setGeneList";
        //    }
        //
        //    public void run() {

        final CursorToken token = WaitCursorManager.showWaitCursor();

        SwingUtilities.invokeLater(new NamedRunnable() {
            public void run() {
                try {
                    if (listID == null) {
                        session.setCurrentGeneList(null);
                    } else {
                        GeneList gl = GeneListManager.getInstance().getGeneList(listID);

                        if (recordHistory) {
                            session.getHistory().push("List: " + listID, 0);
                        }
                        session.setCurrentGeneList(gl);
                    }
                    Preloader.preload();
                    resetFrames();
                } finally {
                    WaitCursorManager.removeWaitCursor(token);

                }
            }

            public String getName() {
                return "Set gene list";
            }
        });
        //  }
        // });


    }

    public void setDefaultFrame(String searchString) {
        FrameManager.setToDefaultFrame(searchString);
        resetFrames();
    }

    public void resetFrames() {
        contentPane.getMainPanel().headerPanelContainer.createHeaderPanels();
        for (TrackPanel tp : getTrackPanels()) {
            tp.createDataPanels();
        }

        contentPane.getCommandBar().setGeneListMode(FrameManager.isGeneListMode());
        contentPane.getMainPanel().revalidate();
        contentPane.getMainPanel().applicationHeaderPanel.revalidate();
        contentPane.getMainPanel().repaint();
    }


    public void enableRemoveGenomes() {

        menuBar.enableRemoveGenomes();

    }


    /**
     * Open the user preferences dialog
     */
    final public void doViewPreferences() {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                boolean originalSingleTrackValue =
                        PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY);

                PreferencesEditor dialog = new PreferencesEditor(mainFrame, true);
                dialog.setVisible(true);


                if (dialog.isCanceled()) {
                    resetStatusMessage();
                    return;

                }


                try {

                    //Should data and feature panels be combined ?
                    boolean singlePanel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY);
                    if (originalSingleTrackValue != singlePanel) {
                        JOptionPane.showMessageDialog(mainFrame, "Panel option change will take affect after restart.");
                    }


                } finally {

                    // Update the state of the current tracks for drawing purposes
                    doRefresh();
                    resetStatusMessage();

                }


            }
        });
    }

    final public void doExitApplication() {

        // Store recent sessions
        if (!getRecentSessionList().isEmpty()) {

            int size = getRecentSessionList().size();
            if (size > UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST) {
                size = UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST;
            }

            String recentSessions = "";
            for (int i = 0; i <
                    size; i++) {
                recentSessions += getRecentSessionList().get(i);

                if (i < (size - 1)) {
                    recentSessions += ";";
                }

            }
            PreferenceManager.getInstance().remove(PreferenceManager.RECENT_SESSION_KEY);
            PreferenceManager.getInstance().setRecentSessions(recentSessions);
        }


    }

    final public void doShowAttributeDisplay(boolean enableAttributeView) {

        boolean oldState = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY);

        // First store the newly requested state
        PreferenceManager.getInstance().setShowAttributeView(enableAttributeView);

        //menuItem.setSelected(enableAttributeView);

        // Now, if the state has actually change we
        // need to refresh everything
        if (oldState != enableAttributeView) {
            doRefresh();
        }


    }


    final public void doRefresh() {

        contentPane.getMainPanel().revalidate();
        mainFrame.repaint();
        //getContentPane().repaint();
    }

    final public void refreshCommandBar() {
        contentPane.getCommandBar().updateCurrentCoordinates();
    }


// TODO -- move all of this attribute stuf out of IGV,  perhaps to

    // some Attribute helper class.

    final public void doSelectDisplayableAttribute() {

        List<String> allAttributes = AttributeManager.getInstance().getAttributeNames();
        Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        final CheckListDialog dlg = new CheckListDialog(mainFrame, allAttributes, hiddenAttributes, false);
        dlg.setVisible(true);

        if (!dlg.isCanceled()) {
            IGV.getInstance().getSession().setHiddenAttributes(dlg.getNonSelections());
            doRefresh();
        }
    }


    final public void saveImage(Component target) {
        saveImage(target, "igv_snapshot");
    }

    final public void saveImage(Component target, String title) {
        contentPane.getStatusBar().setMessage("Creating image...");
        File defaultFile = new File(title + ".png");
        try {
            //createSnapshot(this, defaultFile);
            createSnapshot(target, defaultFile);
        } catch (Exception e) {
            log.error("Error exporting  image ", e);
            MessageUtils.showMessage(("Error encountered while exporting image: " + e.getMessage()));

        } finally {
            resetStatusMessage();

        }
    }

    public boolean isExportingSnapshot() {
        return isExportingSnapshot;
    }

    final public void createSnapshot(final Component target, final File defaultFile) {

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            contentPane.getStatusBar().setMessage("Exporting image: " + defaultFile.getAbsolutePath());
            File file = selectSnapshotFile(defaultFile);
            if (file == null) {
                return;
            }
            createSnapshotNonInteractive(target, file);
        } catch (Exception e) {
            log.error("Error creating exporting image ", e);
            MessageUtils.showMessage(("Error creating the image file: " + defaultFile + "<br> "
                    + e.getMessage()));
        } finally {
            WaitCursorManager.removeWaitCursor(token);
            resetStatusMessage();
        }

    }


    public void createSnapshotNonInteractive(File file) {
        createSnapshotNonInteractive(contentPane.getMainPanel(), file);
    }

    protected void createSnapshotNonInteractive(Component target, File file) {

        log.debug("Creating snapshot: " + file.getName());

        String extension = SnapshotUtilities.getFileExtension(file.getAbsolutePath());

        SnapshotFileType type = SnapshotUtilities.getSnapshotFileType(extension);

        // If valid extension
        if (type != SnapshotFileType.NULL) {

            boolean doubleBuffered = RepaintManager.currentManager(contentPane).isDoubleBufferingEnabled();
            try {
                setExportingSnapshot(true);
                doComponentSnapshot(target, file, type);

            } finally {
                setExportingSnapshot(false);
            }
        }

        log.debug("Finished creating snapshot: " + file.getName());
    }

    public File selectSnapshotFile(File defaultFile) {

        File snapshotDirectory = PreferenceManager.getInstance().getLastSnapshotDirectory();

        JFileChooser fc = new SnapshotFileChooser(snapshotDirectory, defaultFile);
        fc.showSaveDialog(mainFrame);
        File file = fc.getSelectedFile();

        // If a file selection was made
        if (file != null) {
            File directory = file.getParentFile();
            if (directory != null) {
                PreferenceManager.getInstance().setLastSnapshotDirectory(directory);
            }

        }

        return file;
    }


    private void createZoomCursors() throws HeadlessException, IndexOutOfBoundsException {
        if (zoomInCursor == null || zoomOutCursor == null) {
            final Image zoomInImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_IN).getImage();
            final Image zoomOutImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_OUT).getImage();
            final Point hotspot = new Point(10, 10);
            zoomInCursor = mainFrame.getToolkit().createCustomCursor(zoomInImage, hotspot, "Zoom in");
            zoomOutCursor = mainFrame.getToolkit().createCustomCursor(zoomOutImage, hotspot, "Zoom out");

        }

    }

    private void createHandCursor() throws HeadlessException, IndexOutOfBoundsException {
        /*if (handCursor == null) {
            BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

            // Make backgroun transparent
            Graphics2D g = handImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
            g.fill(rect);

            // Draw hand image in middle
            g = handImage.createGraphics();
            g.drawImage(IconFactory.getInstance().getIcon(IconFactory.IconID.OPEN_HAND).getImage(), 0, 0, null);
            handCursor = getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
        }*/

        if (fistCursor == null) {
            BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

            // Make backgroun transparent
            Graphics2D g = handImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(
                    AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
            g.fill(rect);

            // Draw hand image in middle
            g = handImage.createGraphics();
            g.drawImage(IconFactory.getInstance().getIcon(IconFactory.IconID.FIST).getImage(), 0, 0, null);
            fistCursor = mainFrame.getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
        }

    }

    private void createDragAndDropCursor()
            throws HeadlessException, IndexOutOfBoundsException {

        if (dragNDropCursor == null) {
            ImageIcon icon =
                    IconFactory.getInstance().getIcon(
                            IconFactory.IconID.DRAG_AND_DROP);

            int width = icon.getIconWidth();
            int height = icon.getIconHeight();

            BufferedImage dragNDropImage =
                    new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

            // Make background transparent
            Graphics2D g = dragNDropImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(
                    AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, width, height);
            g.fill(rect);

            // Draw DND image
            g =
                    dragNDropImage.createGraphics();
            Image image = icon.getImage();
            g.drawImage(image, 0, 0, null);
            dragNDropCursor =
                    mainFrame.getToolkit().createCustomCursor(
                            dragNDropImage, new Point(0, 0), "Drag and Drop");
        }

    }

    public void resetSession(String sessionPath) {

        LRUCache.clearCaches();

        AttributeManager.getInstance().clearAllAttributes();

        mainFrame.setTitle(UIConstants.APPLICATION_NAME);

        menuBar.resetSessionActions();

        AttributeManager.getInstance().clearAllAttributes();

        SeekableFileStream.closeAllInstances();

        if (session == null) {
            session = new Session(sessionPath);
        } else {
            session.reset(sessionPath);
        }

        contentPane.getMainPanel().resetPanels();

        doRefresh();

    }

    /**
     * Set the status bar message.  If the message equals "Done." intercept
     * and reset to the default "quite" message,  currently the number of tracks
     * loaded.
     *
     * @param message
     */
    public void setStatusBarMessage(String message) {
        if (message.equals("Done.")) {
            resetStatusMessage();
        }

        contentPane.getStatusBar().setMessage(message);
    }

    /**
     * Resets factory settings. this is not the same as reset user defaults
     * DO NOT DELETE used when debugging
     */
    public void resetToFactorySettings() {

        try {
            PreferenceManager.getInstance().clear();
            boolean isShow = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY);
            doShowAttributeDisplay(isShow);
            doRefresh();

        } catch (Exception e) {
            String message = "Failure while resetting preferences!";
            log.error(message, e);
            MessageUtils.showMessage(message + ": " + e.getMessage());
        }

    }

    public void setFilterMatchAll(boolean value) {
        menuBar.setFilterMatchAll(value);
    }

    public boolean isFilterMatchAll() {
        return menuBar.isFilterMatchAll();
    }

    public void setFilterShowAllTracks(boolean value) {
        menuBar.setFilterShowAllTracks(value);

    }

    public boolean isFilterShowAllTracks() {
        return menuBar.isFilterShowAllTracks();
    }

    /**
     * Add a new data panel set
     */
    public TrackPanelScrollPane addDataPanel(String name) {
        return contentPane.getMainPanel().addDataPanel(name);
    }


    /**
     * Return the panel with the given name.  This is called infrequently, and doesn't need to be fast (linear
     * search is fine).
     *
     * @param name
     * @return
     */
    public TrackPanel getTrackPanel(String name) {
        for (TrackPanel sp : getTrackPanels()) {
            if (name.equals(sp.getName())) {
                return sp;
            }
        }

        // If we get this far this is a new panel
        TrackPanelScrollPane sp = addDataPanel(name);
        return sp.getTrackPanel();
    }


    /**
     * Return an ordered list of track panels.  This method is provided primarily for storing sessions, where
     * the track panels need to be stored in order.
     */
    public List<TrackPanel> getTrackPanels() {
        return contentPane.getMainPanel().getTrackPanels();
    }


    public boolean scrollToTrack(String trackName) {
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.getScrollPane().getNamePanel().scrollTo(trackName)) {
                return true;
            }
        }
        return false;
    }


    public Session getSession() {
        return session;
    }

    final public void doRestoreSession(final File sessionFile,
                                       final String locus) {

        if (sessionFile.exists()) {

            doRestoreSession(sessionFile.getAbsolutePath(), locus, false);

        } else {
            String message = "Session file does not exist! : " + sessionFile.getAbsolutePath();
            log.error(message);
            MessageUtils.showMessage(message);
        }

    }


    final public void doRestoreSession(final URL sessionURL,
                                       final String locus) throws Exception {

        doRestoreSession(URLDecoder.decode(sessionURL.toExternalForm(), "UTF-8"), locus, false);

    }

    final public void doRestoreSession(final String sessionPath,
                                       final String locus,
                                       final boolean merge) {

        SwingWorker worker = new SwingWorker() {

            CursorToken cursorToken = null;

            @Override
            protected Object doInBackground() throws Exception {
                InputStream inputStream = null;
                try {

                    cursorToken = WaitCursorManager.showWaitCursor();

                    if (!merge) {
                        // Do this first, it closes all open SeekableFileStreams.
                        resetSession(sessionPath);
                    }

                    setStatusBarMessage("Opening session...");
                    inputStream = ParsingUtils.openInputStream(new ResourceLocator(sessionPath));

                    final SessionReader sessionReader = new SessionReader(IGV.this);

                    sessionReader.loadSession(inputStream, session, sessionPath);

                    String searchText = locus == null ? session.getLocus() : locus;

                    // NOTE: Nothing to do if chr == all
                    if (!FrameManager.isGeneListMode() && searchText != null &&
                            !searchText.equals(Globals.CHR_ALL) && searchText.trim().length() > 0) {
                        goToLocus(searchText);
                    }


                    mainFrame.setTitle(UIConstants.APPLICATION_NAME + " - Session: " + sessionPath);
                    LRUCache.clearCaches();


                    double[] dividerFractions = session.getDividerFractions();
                    if (dividerFractions != null) {
                        contentPane.getMainPanel().setDividerFractions(dividerFractions);
                    }
                    session.clearDividerLocations();

                    //If there's a RegionNavigatorDialog, kill it.
                    //this could be done through the Observer that RND uses, I suppose.  Not sure that's cleaner
                    RegionNavigatorDialog.destroyActiveInstance();

                    if (!getRecentSessionList().contains(sessionPath)) {
                        getRecentSessionList().addFirst(sessionPath);
                    }


                    doRefresh();
                } catch (Exception e) {
                    String message = "Failed to load session! : " + sessionPath;
                    log.error(message, e);
                    MessageUtils.showMessage(message + ": " + e.getMessage());

                } finally {
                    if (inputStream != null) {
                        try {
                            inputStream.close();
                        } catch (IOException iOException) {
                            log.error("Error closing session stream", iOException);
                        }
                        resetStatusMessage();
                    }
                }
                return null;
            }

            @Override
            protected void done() {
                if (cursorToken != null) {
                    WaitCursorManager.removeWaitCursor(cursorToken);
                }
            }
        };

        worker.execute();


    }

    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        contentPane.getStatusBar().setMessage("" +
                getVisibleTrackCount() + " tracks loaded");

    }


    public void rebuildGenomeDropdownList(Set excludedArchivesUrls) {
        contentPane.getCommandBar().rebuildGenomeItemList(excludedArchivesUrls);
    }

    public void showLoadedTrackCount() {
        contentPane.getStatusBar().setMessage("" +
                getVisibleTrackCount() +
                " track(s) currently loaded");
    }

    private void closeWindow(final ProgressBar progressBar) {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                progressBar.close();
            }
        });
    }

    /**
     * Method provided to jump to a locus synchronously.  Used for port command options
     *
     * @param locus
     */
    public void goToLocus(String locus) {

        contentPane.getCommandBar().searchByLocus(locus, false);
    }

    public void tweakPanelDivider() {
        contentPane.getMainPanel().tweakPanelDivider();
    }

    public void removeDataPanel(String name) {
        contentPane.getMainPanel().removeDataPanel(name);
    }

    public void layoutMainPanel() {
        contentPane.getMainPanel().doLayout();
    }

    public MainPanel getMainPanel() {
        return contentPane.getMainPanel();
    }

    public void setExportingSnapshot(boolean exportingSnapshot) {
        isExportingSnapshot = exportingSnapshot;
        if (isExportingSnapshot) {
            RepaintManager.currentManager(contentPane).setDoubleBufferingEnabled(false);
        } else {
            RepaintManager.currentManager(contentPane).setDoubleBufferingEnabled(true);
        }
    }

    public LinkedList<String> getRecentSessionList() {
        return recentSessionList;
    }

    public void setRecentSessionList(LinkedList<String> recentSessionList) {
        this.recentSessionList = recentSessionList;
    }

    public IGVContentPane getContentPane() {
        return contentPane;
    }

    public GenomeManager getGenomeManager() {
        return genomeManager;
    }

    JCheckBoxMenuItem showPeakMenuItem;
    PeakCommandBar peakCommandBar;

    public void addCommandBar(PeakCommandBar cb) {
        this.peakCommandBar = cb;
        contentPane.add(peakCommandBar);
        contentPane.invalidate();

        showPeakMenuItem = new JCheckBoxMenuItem("Show peaks toolbar");
        showPeakMenuItem.setSelected(true);
        showPeakMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                if (showPeakMenuItem.isSelected()) {
                    contentPane.add(peakCommandBar);
                } else {
                    contentPane.remove(peakCommandBar);
                }
            }
        });

        menuBar.getViewMenu().addSeparator();
        menuBar.getViewMenu().add(showPeakMenuItem);
    }

    public boolean isSuppressTooltip() {

        return contentPane != null && contentPane.getCommandBar().isSuppressTooltip();
    }

    public void openStatusWindow() {
        if (statusWindow == null) {
            statusWindow = new StatusWindow();
        }
        statusWindow.setVisible(true);
    }

    public void setStatusWindowText(String text) {
        if (statusWindow != null && statusWindow.isVisible()) {
            statusWindow.updateText(text);
        }
    }


    public void loadResources(Collection<ResourceLocator> locators) {

        //Set<TrackPanel> changedPanels = new HashSet();

        log.info("Loading" + locators.size() + " resources.");
        final MessageCollection messages = new MessageCollection();


        // Load files concurrently -- TODO, put a limit on # of threads?
        List<Thread> threads = new ArrayList<Thread>(locators.size());

        for (final ResourceLocator locator : locators) {

            // If its a local file, check explicitly for existence (rather than rely on exception)
            if (locator.isLocal()) {
                File trackSetFile = new File(locator.getPath());
                if (!trackSetFile.exists()) {
                    messages.append("File not found: " + locator.getPath() + "\n");
                    continue;
                }
            }

            Runnable runnable = new Runnable() {
                public void run() {
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
                                if (t instanceof VariantTrack && ((VariantTrack) t).getAllSamples().size() > 10) {
                                    String newPanelName = "Panel" + System.currentTimeMillis();
                                    panel = addDataPanel(newPanelName).getTrackPanel();
                                }
                            }
                            panel.addTracks(tracks);
                        }
                    } catch (Throwable e) {
                        log.error("Error loading tracks", e);
                        messages.append(e.getMessage());
                    }
                }
            };

            //Thread thread = new Thread(runnable);
            //thread.start();
            //threads.add(thread);
            runnable.run();
        }

        // Wait for all threads to complete
        for (Thread t : threads) {
            try {
                t.join();
            } catch (InterruptedException ignore) {
            }
        }

        resetGroups();
        resetOverlayTracks();

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

        TrackLoader loader = new TrackLoader();
        try {
            List<Track> newTracks = loader.load(locator, this);
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
        doRefresh();
    }

    public Set<TrackType> getLoadedTypes() {
        Set<TrackType> types = new HashSet();
        for (Track t : getAllTracks(false)) {
            TrackType type = t.getTrackType();
            if (t != null) {
                types.add(type);
            }
        }
        return types;
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
            return getTrackPanel(DATA_PANEL_NAME);
        } else if (path.endsWith(".sam") || path.endsWith(".bam") || path.endsWith(".bam.hg19") ||
                path.endsWith(".sam.list") || path.endsWith(".bam.list") ||
                path.endsWith(".aligned") || path.endsWith(".sorted.txt")) {

            String newPanelName = "Panel" + System.currentTimeMillis();
            return addDataPanel(newPanelName).getTrackPanel();
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
            return getTrackPanel(FEATURE_PANEL_NAME);
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
            return getTrackPanel(FEATURE_PANEL_NAME);
        } else {
            return getTrackPanel(DATA_PANEL_NAME);
        }
    }

    public Track getGeneTrack() {
        return geneTrack;
    }


    public SequenceTrack getSequenceTrack() {
        return sequenceTrack;
    }


    public void reset() {
        groupByAttribute = null;
        for (TrackPanel sp : getTrackPanels()) {
            if (DATA_PANEL_NAME.equals(sp.getName())) {
                sp.reset();
                break;
            }
        }
        groupListeners.clear();
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


    /**
     * Reset the overlay tracks collection.  Currently the only overlayable track
     * type is Mutation.  This method finds all mutation tracks and builds a map
     * of key -> mutatinon track,  where the key is the specified attribute value
     * for linking tracks for overlay.
     */
    public void resetOverlayTracks() {
        overlayTracksMap.clear();
        overlaidTracks.clear();


        // Old option to allow overlaying based on an arbitrary attribute.
        // String overlayAttribute = igv.getSession().getOverlayAttribute();

        for (Track track : getAllTracks(false)) {
            if (track != null && track.getTrackType() == TrackType.MUTATION) {

                String sample = track.getSample();

                if (sample != null) {
                    List<Track> trackList = overlayTracksMap.get(sample);

                    if (trackList == null) {
                        trackList = new ArrayList();
                        overlayTracksMap.put(sample, trackList);
                    }

                    trackList.add(track);
                }
            }

        }

        for (Track track : getAllTracks(false)) {
            if (track != null) {  // <= this should not be neccessary
                if (track.getTrackType() != TrackType.MUTATION) {
                    String sample = track.getSample();
                    if (sample != null) {
                        List<Track> trackList = overlayTracksMap.get(sample);
                        if (trackList != null) overlaidTracks.addAll(trackList);
                    }
                }
            }
        }

        boolean displayOverlays = getSession().getOverlayMutationTracks();
        for (Track track : getAllTracks(false)) {
            if (track != null) {
                if (track.getTrackType() == TrackType.MUTATION) {
                    track.setOverlayed(displayOverlays && overlaidTracks.contains(track));
                }
            }
        }
    }


    /**
     * Return tracks overlaid on "track"
     * // TODO -- why aren't overlaid tracks stored in a track member?  This seems unnecessarily complex
     *
     * @param track
     * @return
     */
    public List<Track> getOverlayTracks(Track track) {
        String sample = track.getSample();
        if (sample != null) {
            return overlayTracksMap.get(sample);
        }
        return null;
    }

    public int getVisibleTrackCount() {
        int count = 0;
        for (TrackPanel tsv : getTrackPanels()) {
            count += tsv.getVisibleTrackCount();

        }
        return count;
    }

    /**
     * Return the list of all tracks in the order they appear on the screen
     *
     * @param includeGeneTrack if false exclude gene and reference sequence tracks.
     * @return
     */
    public List<Track> getAllTracks(boolean includeGeneTrack) {
        List<Track> allTracks = new ArrayList<Track>();

        for (TrackPanel tp : getTrackPanels()) {
            allTracks.addAll(tp.getTracks());
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
     * Return the complete set of unique DataResourceLocators currently loaded
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


    public void setAllTrackHeights(int newHeight) {
        for (Track track : getAllTracks(false)) {
            track.setHeight(newHeight);
        }

    }


    public void removeTracks(Collection<Track> tracksToRemove) {

        // Make copy of list as we will be modifying the original in the loop
        List<TrackPanel> panels = getTrackPanels();
        for (TrackPanel trackPanel : panels) {
            trackPanel.removeTracks(tracksToRemove);

            if (!trackPanel.hasTracks()) {
                removeDataPanel(trackPanel.getName());
            }
        }

        for (Track t : tracksToRemove) {
            if (t instanceof DragListener) {
                DragEventManager.getInstance().removeDragListener((DragListener) t);
            }
            if (t instanceof TrackGroupEventListener) {
                removeGroupEventListener((TrackGroupEventListener) t);
            }
            if(t instanceof AlignmentTrackEventListener) {
                removeAlignmentTrackEvent((AlignmentTrackEventListener) t);
            }
        }
    }


    /**
     * @param reader        a reader for the gene (annotation) file.
     * @param genome
     * @param geneFileName
     * @param geneTrackName
     */
    public void createGeneTrack(Genome genome, AsciiLineReader reader, String geneFileName, String geneTrackName,
                                String annotationURL) {

        FeatureDB.clearFeatures();
        FeatureTrack geneFeatureTrack = null;

        if (reader != null) {
            FeatureParser parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(geneFileName), genome);
            if (parser == null) {
                MessageUtils.showMessage("ERROR: Unrecognized annotation file format: " + geneFileName +
                        "<br>Annotations for genome: " + genome.getId() + " will not be loaded.");
            } else {
                List<org.broad.tribble.Feature> genes = parser.loadFeatures(reader);
                String name = geneTrackName;
                if (name == null) name = "Genes";

                String id = genome.getId() + "_genes";
                geneFeatureTrack = new FeatureTrack(id, name, new FeatureCollectionSource(genes, genome));
                geneFeatureTrack.setMinimumHeight(5);
                geneFeatureTrack.setHeight(35);
                geneFeatureTrack.setPreferredHeight(35);
                geneFeatureTrack.setRendererClass(IGVFeatureRenderer.class);
                geneFeatureTrack.setColor(Color.BLUE.darker());
                TrackProperties props = parser.getTrackProperties();
                if (props != null) {
                    geneFeatureTrack.setProperties(parser.getTrackProperties());
                }
                geneFeatureTrack.setUrl(annotationURL);
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
        for (TrackPanel tsv : getTrackPanels()) {
            foundSeqTrack = tsv.replaceTrack(sequenceTrack, newSeqTrack);
            if (foundSeqTrack) {
                break;
            }
        }

        boolean foundGeneTrack = false;
        for (TrackPanel tsv : getTrackPanels()) {
            foundGeneTrack = tsv.replaceTrack(geneTrack, newGeneTrack);
            if (foundGeneTrack) {
                break;
            }
        }


        if (!foundGeneTrack || !foundSeqTrack) {
            TrackPanel panel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY) ?
                    getTrackPanel(DATA_PANEL_NAME) : getTrackPanel(FEATURE_PANEL_NAME);

            if (!foundSeqTrack) panel.addTrack(newSeqTrack);
            if (!foundGeneTrack && newGeneTrack != null) panel.addTrack(newGeneTrack);

        }

        // Keep a reference to this track so it can be removed
        geneTrack = newGeneTrack;
        sequenceTrack = newSeqTrack;

    }


    /////////////////////////////////////////////////////////////////////////////////////////
    // Sorting


    /**
     * Sort all groups (data and feature) by attribute value(s).  Tracks are
     * sorted within groups.
     *
     * @param attributeNames
     * @param ascending
     */
    public void sortAllTracksByAttributes(final String attributeNames[], final boolean[] ascending) {
        assert attributeNames.length == ascending.length;

        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.sortTracksByAttributes(attributeNames, ascending);
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

        // Create a rank order of samples.  This is done globally so sorting is consistent across groups and panels.
        final List<String> sortedSamples = sortSamplesByRegionScore(region, type, frame);

        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.sortByRegionsScore(r, type, frame, sortedSamples);
        }
        repaintDataPanels();
    }


    /**
     * Sort a collection of tracks by a score over a region.
     *
     * @param region
     * @param type
     * @param frame
     */
    private List<String> sortSamplesByRegionScore(final RegionOfInterest region,
                                                  final RegionScoreType type,
                                                  final ReferenceFrame frame) {

        // Get the sortable tracks for this score (data) type
        final List<Track> allTracks = getAllTracks(false);
        final List<Track> tracksWithScore = new ArrayList(allTracks.size());
        for (Track t : allTracks) {
            if (t.isRegionScoreType(type)) {
                tracksWithScore.add(t);
            }
        }

        // Sort the "sortable" tracks
        sortByRegionScore(tracksWithScore, region, type, frame);

        // Now get sample order from sorted tracks, use to sort (tracks which do not implement the selected "sort by" score)
        List<String> sortedSamples = new ArrayList(tracksWithScore.size());
        for (Track t : tracksWithScore) {
            String att = t.getSample(); //t.getAttributeValue(linkingAtt);
            if (att != null) {
                sortedSamples.add(att);
            }

        }

        return sortedSamples;
    }

    private void sortByRegionScore(List<Track> tracks,
                                   final RegionOfInterest region,
                                   final RegionScoreType type,
                                   final ReferenceFrame frame) {
        if ((tracks != null) && (region != null) && !tracks.isEmpty()) {
            final int zoom = Math.max(0, frame.getZoom());
            final String chr = region.getChr();
            final int start = region.getStart();
            final int end = region.getEnd();

            Comparator<Track> c = new Comparator<Track>() {

                public int compare(Track t1, Track t2) {
                    try {
                        if (t1 == null && t2 == null) return 0;
                        if (t1 == null) return 1;
                        if (t2 == null) return -1;

                        float s1 = t1.getRegionScore(chr, start, end, zoom, type, frame);
                        float s2 = t2.getRegionScore(chr, start, end, zoom, type, frame);

                        if (s2 > s1) {
                            return 1;
                        } else if (s1 < s2) {
                            return -1;
                        } else {
                            return 0;
                        }
                    } catch (Exception e) {
                        log.error("Error sorting tracks. Sort might not be accurate.", e);
                        return 0;
                    }

                }
            };
            Collections.sort(tracks, c);

        }
    }


    ////////////////////////////////////////////////////////////////////////////////////////
    // Groups

    public String getGroupByAttribute() {
        return groupByAttribute;
    }


    public void setGroupByAttribute(String attributeName) {
        groupByAttribute = attributeName;
        resetGroups();
        // Some tracks need to respond to changes in grouping, fire notification event
        notifyGroupEvent();
    }


    private void resetGroups() {
        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.groupTracksByAttribute(groupByAttribute);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Events

    public synchronized void addGroupEventListener(TrackGroupEventListener l) {
        groupListeners.add(new SoftReference<TrackGroupEventListener>(l));
    }

    public synchronized void removeGroupEventListener(TrackGroupEventListener l) {

        for (Iterator<SoftReference<TrackGroupEventListener>> it = groupListeners.iterator(); it.hasNext(); ) {
            TrackGroupEventListener listener = it.next().get();
            if (listener != null && listener == l) {
                it.remove();
                break;
            }
        }
    }

    public void notifyGroupEvent() {
        TrackGroupEvent e = new TrackGroupEvent(this);
        for (SoftReference<TrackGroupEventListener> ref : groupListeners) {
            TrackGroupEventListener l = ref.get();
            l.onTrackGroupEvent(e);
        }
    }

    public synchronized void addAlignmentTrackEventListener(AlignmentTrackEventListener l) {
        alignmentTrackListeners.add(new SoftReference<AlignmentTrackEventListener>(l));
    }

    public synchronized void removeAlignmentTrackEvent(AlignmentTrackEventListener l) {
        for (Iterator<SoftReference<AlignmentTrackEventListener>> it = alignmentTrackListeners.iterator(); it.hasNext(); ) {
            AlignmentTrackEventListener listener = it.next().get();
            if (listener != null && listener == l) {
                it.remove();
                break;
            }
        }
    }

    public void notifyAlignmentTrackEvent(Object source, AlignmentTrackEvent.Type type) {
        AlignmentTrackEvent e = new AlignmentTrackEvent(source, type);
        for (SoftReference<AlignmentTrackEventListener> ref : alignmentTrackListeners) {
            AlignmentTrackEventListener l = ref.get();
            l.onAlignmentTrackEvent(e);
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////
    // Startup


    public void startUp(Main.IGVArgs igvArgs) {

        if (log.isDebugEnabled()) {
            log.debug("startUp");
        }

        SwingWorker worker = new StartupWorker(igvArgs);
        worker.execute();
    }

    public boolean isStartupComplete() {
        return startupComplete;
    }

    /**
     * Swing worker class to startup IGV
     */
    public class StartupWorker extends SwingWorker {
        Main.IGVArgs igvArgs;

        StartupWorker(Main.IGVArgs args) {
            this.igvArgs = args;

        }

        /**
         * Do the actual work
         *
         * @return
         * @throws Exception
         */
        @Override
        protected Object doInBackground() throws Exception {


            final ProgressMonitor monitor = new ProgressMonitor();
            final ProgressBar progressBar = ProgressBar.showProgressDialog(mainFrame, "Initializing...", monitor, false);
            monitor.fireProgressChange(20);

            try {
                contentPane.getCommandBar().initializeGenomeList(monitor);
            } catch (FileNotFoundException ex) {
                JOptionPane.showMessageDialog(mainFrame, "Error initializing genome list: " + ex.getMessage());
                log.error("Error initializing genome list: ", ex);
            } catch (NoRouteToHostException ex) {
                JOptionPane.showMessageDialog(mainFrame, "Network error initializing genome list: " + ex.getMessage());
                log.error("Network error initializing genome list: ", ex);
            } finally {
                monitor.fireProgressChange(50);
                closeWindow(progressBar);
            }

            final PreferenceManager preferenceManager = PreferenceManager.getInstance();
            if (igvArgs.getGenomeId() != null) {
                selectGenomeFromList(igvArgs.getGenomeId());
            } else if (igvArgs.getSessionFile() == null) {
                String genomeId = preferenceManager.getDefaultGenome();
                contentPane.getCommandBar().selectGenomeFromList(genomeId);
            }

            //If there is an argument assume it is a session file or url
            if (igvArgs.getSessionFile() != null || igvArgs.getDataFileString() != null) {

                if (log.isDebugEnabled()) {
                    log.debug("Loadding session data");
                }

                final IndefiniteProgressMonitor indefMonitor = new IndefiniteProgressMonitor(60);
                final ProgressBar bar2 = ProgressBar.showProgressDialog(mainFrame, "Loading session data", indefMonitor, false);
                indefMonitor.start();

                try {

                    if (log.isDebugEnabled()) {
                        log.debug("Calling restore session");
                    }


                    if (igvArgs.getSessionFile() != null) {
                        if (HttpUtils.getInstance().isURL(igvArgs.getSessionFile())) {
                            URL url = new URL(igvArgs.getSessionFile());
                            doRestoreSession(url, igvArgs.getLocusString());
                        } else {
                            File sf = new File(igvArgs.getSessionFile());
                            if (sf.exists()) {
                                doRestoreSession(sf, igvArgs.getLocusString());
                            }
                        }
                    } else if (igvArgs.getDataFileString() != null) {
                        // Not an xml file, assume its a list of data files
                        String[] tokens = igvArgs.getDataFileString().split(",");
                        String indexFile = igvArgs.getIndexFile();
                        List<ResourceLocator> locators = new ArrayList();
                        for (String p : tokens) {
                            if (FileUtils.isRemote(p)) {
                                p = URLDecoder.decode(p);
                            }
                            ResourceLocator rl = new ResourceLocator(p);
                            if (igvArgs.getName() != null) {
                                rl.setName(igvArgs.getName());
                            }
                            rl.setIndexPath(indexFile);
                            locators.add(rl);
                        }
                        loadResources(locators);
                    }


                } catch (Exception ex) {
                    String tmp = igvArgs.getSessionFile() != null ? igvArgs.getSessionFile() : igvArgs.getDataFileString();
                    JOptionPane.showMessageDialog(mainFrame, "<html>Error loading session: " + tmp + "<br>" + ex.toString());
                    log.error("Error loading session: " + tmp, ex);

                    // Session load failed, load default genome
                    String genomeId = preferenceManager.getDefaultGenome();
                    contentPane.getCommandBar().selectGenomeFromList(genomeId);

                }


                indefMonitor.stop();
                closeWindow(bar2);
            }

            if (igvArgs.getLocusString() != null) {
                goToLocus(igvArgs.getLocusString());
            }

            session.recordHistory();

            // Start up a port listener.  Port # can be overriden with "-p" command line switch
            boolean portEnabled = preferenceManager.getAsBoolean(PreferenceManager.PORT_ENABLED);
            String portString = igvArgs.getPort();
            if (portEnabled || portString != null) {
                // Command listner thread
                int port = preferenceManager.getAsInt(PreferenceManager.PORT_NUMBER);
                if (portString != null) {
                    port = Integer.parseInt(portString);
                }
                CommandListener.start(port);
            }


            startupComplete = true;

            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    mainFrame.setVisible(true);
                }
            });


            return null;
        }


        /**
         * Called when the background thread is complete (IGV window is open and data loaded).
         */
        @Override
        protected void done() {
            if (igvArgs.getBatchFile() != null) {
                LongRunningTask.submit(new BatchRunner(igvArgs.getBatchFile()));
            }
        }
    }


}
