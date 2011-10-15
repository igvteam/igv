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
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.*;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.lists.Preloader;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.peaks.PeakCommandBar;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionReader;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.TrackManager;

import static org.broad.igv.ui.WaitCursorManager.CursorToken;

import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.filefilters.CoverageFileFilter;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;

import static org.broad.igv.ui.util.SnapshotUtilities.*;

import org.broad.igv.ui.util.ProgressMonitor;

import static org.broad.igv.ui.util.UIUtilities.getFileChooser;

import org.broad.igv.ui.filefilters.AlignmentFileFilter;

import org.broad.igv.util.*;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.*;
import java.net.*;
import java.util.*;
import java.util.List;
import java.util.prefs.Preferences;

/**
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

    // Manager classes
    private TrackManager trackManager;
    private GenomeManager genomeManager;

    // FileChooser Dialogs
    private FileChooserDialog trackFileChooser;
    private FileChooser snapshotFileChooser;


    // Misc state
    private LinkedList<String> recentSessionList = new LinkedList<String>();
    private boolean isExportingSnapshot = false;
    private boolean startupComplete = false;


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
                ToolTipManager.sharedInstance().setEnabled(false);
                IGVPopupMenu.closeAll();
            }


            @Override
            public void windowDeactivated(WindowEvent windowEvent) {
                ToolTipManager.sharedInstance().setEnabled(false);
                IGVPopupMenu.closeAll();
            }

            @Override
            public void windowActivated(WindowEvent windowEvent) {
                ToolTipManager.sharedInstance().setEnabled(true);
            }

            @Override
            public void windowGainedFocus(WindowEvent windowEvent) {
                ToolTipManager.sharedInstance().setEnabled(true);
            }
        });


        session = new Session(null);
        trackManager = new TrackManager(this);

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
        contentPane = new IGVContentPane(trackManager);
        menuBar = new IGVMenuBar();

        rootPane.setContentPane(contentPane);
        rootPane.setJMenuBar(menuBar);
        glassPane = rootPane.getGlassPane();
        glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        glassPane.addMouseListener(new MouseAdapter() {
        });
        dNdGlassPane = new GhostGlassPane();

        mainFrame.pack();

        // TODO -- refactor to eliminate these
        initializeSnapshot();
        initializeDialogs();

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

    public void setSelectedRegion(RegionOfInterest region) {
        //if (region != regionOfInterestPane.getSelectedRegion()) {
        //    regionOfInterestPane.setSelectedRegion(region);
        //    repaintDataPanels();
        //}
    }


    public Dimension getPreferredSize() {
        return UIConstants.preferredSize;
    }


    // TODO -- eliminate this shared file chooser,  and all "shared" dialogs like this.

    private void initializeDialogs() {

        // Create Track Chooser
        //  Note --  why are these reused ? (JTR)
        trackFileChooser = new FileChooserDialog(mainFrame, true);
        trackFileChooser.addChoosableFileFilter(new AlignmentFileFilter());
        trackFileChooser.addChoosableFileFilter(new CoverageFileFilter());

        // This hack is ugly, but I can't see any other way to set the default file filter to "All"
        trackFileChooser.setFileFilter(trackFileChooser.getChoosableFileFilters()[0]);


    }

    public FileChooserDialog getTrackFileChooser() {
        return trackFileChooser;
    }

    private void initializeSnapshot() {

        File snapshotDirectory = PreferenceManager.getInstance().getLastSnapshotDirectory();


        // File Filters
        FileFilter[] fileFilters = SnapshotUtilities.getAllSnapshotFileFilters();

        snapshotFileChooser = getFileChooser(snapshotDirectory, null, fileFilters);
        snapshotFileChooser.setDialogTitle("Snapshot File");

        snapshotFileChooser.addPropertyChangeListener(
                new PropertyChangeListener() {

                    public void propertyChange(PropertyChangeEvent e) {

                        File oldFile = null;
                        String property = e.getPropertyName();
                        if (JFileChooser.SELECTED_FILE_CHANGED_PROPERTY.equals(property)) {
                            oldFile = (File) e.getOldValue();
                            snapshotFileChooser.setPreviousFile(oldFile);
                        } else if (JFileChooser.FILE_FILTER_CHANGED_PROPERTY.equals(property)) {

                            if (e.getOldValue() instanceof SnapshotFileFilter &&
                                    e.getNewValue() instanceof SnapshotFileFilter) {

                                SnapshotFileFilter newFilter =
                                        (SnapshotFileFilter) e.getNewValue();

                                File currentDirectory = snapshotFileChooser.getCurrentDirectory();
                                File previousFile = snapshotFileChooser.getPreviousFile();
                                if (previousFile != null) {

                                    File file = null;
                                    if (currentDirectory != null) {
                                        file = new File(currentDirectory, previousFile.getName());
                                    } else {
                                        file = previousFile;
                                    }

                                    final File selectedFile = Utilities.changeFileExtension(
                                            file, newFilter.getExtension());

                                    UIUtilities.invokeOnEventThread(new Runnable() {

                                        public void run() {
                                            snapshotFileChooser.setSelectedFile(selectedFile);
                                            snapshotFileChooser.validate();
                                        }
                                    });
                                }

                            }
                        }
                    }
                });
    }

    public void addRegionOfInterest(RegionOfInterest roi) {
        session.addRegionOfInterestWithNoListeners(roi);
        RegionOfInterestPanel.setSelectedRegion(roi);
        doRefresh();
    }

    void beginROI(JButton button) {
        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
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

        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
            DataPanelContainer dp = tsv.getDataPanel();
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
        if (Globals.batch) {
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
     *
     * @deprecated
     */
    public void repaintDataPanels() {
        repaintDataAndHeaderPanels(false);
    }

    public void repaintNamePanels() {
        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
            tsv.getNamePanel().repaint();
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

        CursorToken token = null;

        try {
            token = WaitCursorManager.showWaitCursor();
            if (locators != null && !locators.isEmpty()) {

                // get current track count per panel.  Needed to detect which panels
                // changed.  Also record panel sizes
                final HashMap<TrackPanelScrollPane, Integer> trackCountMap = new HashMap();
                final HashMap<TrackPanelScrollPane, Integer> panelSizeMap = new HashMap();
                final Collection<TrackPanelScrollPane> scrollPanes = trackManager.getTrackPanelScrollPanes();
                for (TrackPanelScrollPane sp : scrollPanes) {
                    trackCountMap.put(sp, sp.getDataPanel().getAllTracks().size());
                    panelSizeMap.put(sp, sp.getDataPanel().getHeight());
                }

                getTrackManager().loadResources(locators);

                double totalHeight = 0;
                for (TrackPanelScrollPane sp : scrollPanes) {
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


                    }
                });
            }
        } catch (Exception e) {
            if (!(e instanceof ConcurrentModificationException)) {
                if (e.getMessage() != null && e.getMessage().length() > 8) {
                    MessageUtils.showMessage(e.getMessage());
                } else {
                    log.error(e);
                    MessageUtils.showMessage("An error occurred while loading tracks. " +
                            "Please check the logs for details.");
                }
            }
        } finally {
            showLoadedTrackCount();
            if (token != null) {
                WaitCursorManager.removeWaitCursor(token);
            }
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
        for (TrackPanelScrollPane tp : trackManager.getTrackPanelScrollPanes()) {
            tp.getTrackPanel().createDataPanels();
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
            HttpUtils.getInstance().shutdown();
            CommandListener.halt();
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

        // Use default extension if file has none
        if (extension == null) {

            FileFilter filter = snapshotFileChooser.getFileFilter();

            // Figure out the proper extension
            if (!(filter instanceof SnapshotFileFilter)) {
                extension = SnapshotFileType.PNG.getExtension();
            } else {
                extension = ((SnapshotFileFilter) filter).getExtension();
            }

            file = new File((file.getAbsolutePath() + extension));
        }

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

    public File selectSnapshotFile(
            File defaultFile) {

        SnapshotFileFilter snapshotFileFilter = null;
        if (defaultFile != null) {

            String fileExtension = SnapshotUtilities.getFileExtension(defaultFile.getAbsolutePath());
            snapshotFileFilter = SnapshotUtilities.getSnapshotFileFilterForType(
                    SnapshotUtilities.getSnapshotFileType(fileExtension));
        }

        snapshotFileChooser.setFileFilter(snapshotFileFilter);
        snapshotFileChooser.setSelectedFile(defaultFile);

        // Display the dialog
        snapshotFileChooser.showSaveDialog(mainFrame);

        resetStatusMessage();

        File file = snapshotFileChooser.getSelectedFile();

        // If a file selection was made
        if (file != null) {

            File directory = snapshotFileChooser.getCurrentDirectory();
            if (directory != null) {
                PreferenceManager.getInstance().setLastSnapshotDirectory(
                        directory);
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

    public void createNewSession(String sessionName) {

        LRUCache.clearCaches();

        AttributeManager.getInstance().clearAllAttributes();

        mainFrame.setTitle(UIConstants.APPLICATION_NAME);

        menuBar.resetSessionActions();

        AttributeManager.getInstance().clearAllAttributes();

        session = new Session(sessionName);

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
            MessageUtils.showAndLogErrorMessage(mainFrame, message, log, e);
        }

    }

    public void updateTrackState() {

        doRefresh();
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


    public TrackPanel getDataPanel(String name) {
        TrackPanelScrollPane sp = trackManager.getScrollPane(name);
        if (sp == null) {
            sp = addDataPanel(name);
            trackManager.putScrollPane(name, sp);
        }
        return sp.getTrackPanel();
    }


    public boolean scrollToTrack(String trackName) {
        for (TrackPanelScrollPane sp : trackManager.getTrackPanelScrollPanes()) {
            if (sp.getNamePanel().scrollTo(trackName)) {
                return true;
            }

        }
        return false;
    }

    /**
     * Return an ordered list of track panels.  This method is provided primarily for storing sessions, where
     * the track panels need to be stored in order.
     */
    public List<TrackPanel> getTrackPanels() {
        return contentPane.getMainPanel().getTrackPanels();
    }


    public Session getSession() {
        return session;
    }

    final public void doRestoreSession(final File sessionFile,
                                       final String locus) {

        String filePath = "";
        if (sessionFile != null) {

            log.debug("Run doRestoreSession");

            InputStream inputStream = null;
            CursorToken token = WaitCursorManager.showWaitCursor();
            try {
                inputStream = new BufferedInputStream(new FileInputStream(sessionFile));
                doRestoreSession(inputStream, sessionFile.getAbsolutePath(), locus, false);

                String sessionFilePath = sessionFile.getAbsolutePath();
                if (!getRecentSessionList().contains(sessionFilePath)) {
                    getRecentSessionList().addFirst(sessionFilePath);
                }

            } catch (Exception e) {
                String message = "Failed to load session! : " + sessionFile.getAbsolutePath();
                MessageUtils.showAndLogErrorMessage(mainFrame, message, log, e);
            } finally {
                WaitCursorManager.removeWaitCursor(token);
                if (inputStream != null) {
                    try {
                        inputStream.close();
                    } catch (IOException iOException) {
                        log.error("Error closing session stream", iOException);
                    }
                }
            }
            log.debug("Finish doRestoreSession");


        } else {
            String message = "Session file does not exist! : " + filePath;
            MessageUtils.showAndLogErrorMessage(mainFrame, message, log);
        }

    }


    final public void doRestoreSession(final URL sessionURL,
                                       final String locus) throws Exception {

        log.debug("Enter doRestoreSession: " + sessionURL + " " + locus);

        InputStream inputStream = null;
        try {
            inputStream = HttpUtils.getInstance().openConnectionStream(sessionURL);
            doRestoreSession(inputStream, URLDecoder.decode(sessionURL.toExternalForm(), "UTF-8"), locus, false);
        } catch (Exception e) {
            log.error("Error restoring session", e);
        } finally {

            if (inputStream != null) {
                try {
                    inputStream.close();
                } catch (IOException iOException) {
                    log.error("Error closing session stream", iOException);
                }
            }
        }

        log.debug("Exit doRestoreSession");

    }

    final public void doRestoreSession(final InputStream inputStream,
                                       final String sessionPath,
                                       final String locus,
                                       boolean merge) {

        try {
            setStatusBarMessage("Opening session...");

            if (!merge) {
                createNewSession(sessionPath);
            }

            final SessionReader sessionReader = new SessionReader(this);

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

            doRefresh();
        } finally {

            resetStatusMessage();
        }

    }

    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        contentPane.getStatusBar().setMessage("" +
                IGV.getInstance().getTrackManager().getVisibleTrackCount() + " tracks loaded");

    }


    public void rebuildGenomeDropdownList(Set excludedArchivesUrls) {
        contentPane.getCommandBar().rebuildGenomeItemList(excludedArchivesUrls);
    }

    public void showLoadedTrackCount() {
        contentPane.getStatusBar().setMessage("" +
                IGV.getInstance().getTrackManager().getVisibleTrackCount() +
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


    public TrackManager getTrackManager() {
        return trackManager;
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
                            rl.setIndexPath(indexFile);
                            locators.add(rl);
                        }
                        getTrackManager().loadResources(locators);
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
