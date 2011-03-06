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

/*
 * IGVMainFrame.java
 *
 * Created on September 28, 2007, 2:14 PM
 */
package org.broad.igv.ui;

import com.jidesoft.status.StatusBar;
import com.jidesoft.swing.JideBoxLayout;
import com.jidesoft.swing.JideSplitPane;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.IGVDatasetParser;
import org.broad.igv.event.StatusChangeEvent;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.GenomeBuilderDialog;
import org.broad.igv.feature.genome.GenomeDescriptor;
import org.broad.igv.feature.genome.GenomeManager.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.lists.GeneListManagerUI;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.tools.ui.CoverageGui;
import org.broad.igv.tools.ui.IndexGui;
import org.broad.igv.listener.StatusListener;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionReader;
import org.broad.igv.tools.IgvToolsGui;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.TrackManager;

import static org.broad.igv.ui.UIConstants.*;
import static org.broad.igv.ui.WaitCursorManager.CursorToken;

import org.broad.igv.ui.action.*;
import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.event.GlobalKeyDispatcher;
import org.broad.igv.ui.legend.LegendDialog;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;

import static org.broad.igv.ui.util.SnapshotUtilities.*;

import org.broad.igv.ui.util.ProgressMonitor;

import static org.broad.igv.ui.util.UIUtilities.getFileChooser;

import org.broad.igv.ui.filefilters.AlignmentFileFilter;

import org.broad.igv.util.*;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.plaf.basic.BasicBorders;
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

/**
 * @author jrobinso
 */
public class IGVMainFrame extends javax.swing.JFrame {

    private static Logger log = Logger.getLogger(IGVMainFrame.class);
    private static IGVMainFrame theInstance;

    // Cursors
    public static Cursor handCursor;
    public static Cursor fistCursor;
    public static Cursor zoomInCursor;
    public static Cursor zoomOutCursor;
    public static Cursor dragNDropCursor;

    //Session session;
    Session session;
    /**
     * Helper class for managing tracks
     */
    private TrackManager trackManager;

    // Panels, UI elements
    private JMenu extrasMenu;
    private IGVCommandBar igvCommandBar;
    private MainPanel mainPanel;
    private StatusBar statusBar;
    private RemoveUserDefinedGenomeMenuAction removeImportedGenomeAction;
    //TODO -- A lot of state is passed between the embedded filter in this
    // action and the session during save and restore.  Refactor to pass
    // this state in a single object to/from the session.
    private FilterTracksMenuAction filterTracksAction;

    // FileChooser Dialogs
    private FileChooserDialog trackFileChooser;
    private FileChooser snapshotFileChooser;
    private FileChooser genomeImportFileChooser;


    // Glass panes
    Component glassPane;
    GhostGlassPane dNdGlassPane;


    // Misc state
    private LinkedList<String> recentSessionList = new LinkedList<String>();
    private boolean isExportingSnapshot = false;


    // Tracksets
    //private final Map<String, TrackPanelScrollPane> trackSetScrollPanes = new Hashtable();


    public static IGVMainFrame getInstance() {
        if (theInstance == null) {
            if (Globals.isHeadless()) {
                System.err.println("Attempt to instantiate IGVMainFrame in headless mode");
            } else {
                theInstance = new IGVMainFrame();
            }
        }
        return theInstance;
    }

    public static boolean hasInstance() {
        return theInstance != null;
    }


    /**
     * Creates new form IGVMainFrame
     */
    private IGVMainFrame() {

        theInstance = this;
        session = new Session(null);
        trackManager = new TrackManager(this);

        // Create cursors
        createHandCursor();
        createZoomCursors();
        createDragAndDropCursor();

        // Create components
        setTitle(UIConstants.APPLICATION_NAME);
        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        igvCommandBar = new IGVCommandBar(this);
        igvCommandBar.setMinimumSize(new Dimension(250, 33));
        mainPanel = new MainPanel(trackManager);
        getContentPane().add(mainPanel, java.awt.BorderLayout.CENTER);
        createMenuAndToolbar();

        statusBar = new ApplicationStatusBar();
        statusBar.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        getContentPane().add(statusBar, java.awt.BorderLayout.SOUTH);
        pack();

        // Create glass panes
        glassPane = getGlassPane();
        getGlassPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        getGlassPane().addMouseListener(new MouseAdapter() {
        });
        dNdGlassPane = new GhostGlassPane();

        // Add listeners
        // Closing the window should exit the application
        addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                IGVMainFrame.this.doExitApplication();
            }
        });


        // Application initialization
        Runtime.getRuntime().addShutdownHook(new ShutdownThread());
        // TODO -- get these from user preferences
        ToolTipManager.sharedInstance().setDismissDelay(Integer.MAX_VALUE);
        //ToolTipManager.sharedInstance().setReshowDelay(your time in ms);
        //ToolTipManager.sharedInstance().setInitialDelay(your time in ms);

        // TODO -- refactor to eliminate these
        initializeSnapshot();
        initializeDialogs();

        // Anti alias settings.   TODO = Are these neccessary anymore ?
        System.setProperty("awt.useSystemAAFontSettings", "on");
        System.setProperty("swing.aatext", "true");

        // Set the application's previous location and size
        Rectangle applicationBounds = PreferenceManager.getInstance().getApplicationFrameBounds();
        Dimension screenBounds = Toolkit.getDefaultToolkit().getScreenSize();
        if (applicationBounds != null &&
                applicationBounds.getMaxX() < screenBounds.getWidth() &&
                applicationBounds.getMaxY() < screenBounds.getHeight()) {
            setBounds(applicationBounds);
        }

        KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new GlobalKeyDispatcher());
        IGVHttpUtils.updateProxySettings();



    }

    public GhostGlassPane getDnDGlassPane() {
        return dNdGlassPane;
    }

    public void startDnD() {
        setGlassPane(dNdGlassPane);
        getGlassPane().setVisible(true);
    }

    public void endDnD() {
        setGlassPane(glassPane);
        getGlassPane().setVisible(false);
    }

    public void setSelectedRegion(RegionOfInterest region) {
        //if (region != regionOfInterestPane.getSelectedRegion()) {
        //    regionOfInterestPane.setSelectedRegion(region);
        //    repaintDataPanels();
        //}
    }


    @Override
    public Dimension getPreferredSize() {
        return UIConstants.preferredSize;
    }



    // TODO -- eliminate this shared file chooser,  and all "shared" dialogs like this.

    private void initializeDialogs() {

        // Create Track Chooser
        //  Note --  why are these reused ? (JTR)
        trackFileChooser = new FileChooserDialog(this, true);
        trackFileChooser.addChoosableFileFilter(new AlignmentFileFilter());

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


    public void chromosomeChangeEvent() {
        chromosomeChangeEvent(true);
    }

    public void chromosomeChangeEvent(boolean updateCommandBar) {
        igvCommandBar.chromosomeChanged();
        repaintDataAndHeaderPanels(updateCommandBar);

    }

    /**
     * Repaint panels containing data, specifically the dataTrackPanel,
     * featureTrackPanel, and headerPanel.
     */
    public void repaintDataAndHeaderPanels() {
        repaintDataAndHeaderPanels(true);
    }

    public void repaintDataAndHeaderPanels(boolean updateCommandBar) {
        mainPanel.repaint();
        if (updateCommandBar) {
            igvCommandBar.updateCurrentCoordinates();
        }
    }

    public void repaintDataPanels() {
        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
            tsv.getDataPanel().repaint();
        }

    }

    public void repaintNamePanels() {
        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
            tsv.getNamePanel().repaint();
        }

    }

    public void repaintStatusAndZoomSlider() {
        igvCommandBar.repaint();
    }


    public void selectGenomeFromList(String genome) {
        try {
            igvCommandBar.selectGenomeFromList(genome);
        } catch (FileNotFoundException e) {
            log.error("File not found while intializing genome!", e);
        } catch (NoRouteToHostException e) {
            log.error("Error while intializing genome!", e);
        }

    }


    private void createMenuAndToolbar() {

        // Setup the menus
        setJMenuBar(MenuAndToolbarUtils.createMenuBar(createMenus()));

        // Setup the toolbar panel.
        JPanel toolbarPanel = new JPanel();
        toolbarPanel.setBorder(new BasicBorders.MenuBarBorder(Color.GRAY, Color.GRAY));

        getContentPane().add(toolbarPanel, BorderLayout.NORTH);
        toolbarPanel.setLayout(new JideBoxLayout(toolbarPanel));

        // Nothing for this toolbar yet, basically used as a space
        //JPanel namePanelToolBar = new JPanel();
        //namePanelToolBar.setLayout(new JideBoxLayout(namePanelToolBar));
        //namePanelToolBar.setPreferredSize(new Dimension(180, 10));
        //toolbarPanel.add(namePanelToolBar, JideBoxLayout.FLEXIBLE);
        toolbarPanel.add(igvCommandBar, JideBoxLayout.VARY);
    }


    public void doDefineGenome(ProgressMonitor monitor) {

        ProgressBar bar = null;
        File archiveFile = null;

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            GenomeBuilderDialog genomeBuilderDialog =
                    new GenomeBuilderDialog(this, true);

            genomeBuilderDialog.setVisible(true);
            if (genomeBuilderDialog.isCanceled()) {
                return;
            }

            if (monitor != null) {
                bar = ProgressBar.showProgressDialog(IGVMainFrame.getInstance(),
                        "Defining Genome...", monitor, false);
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

            GenomeListItem genomeListItem = GenomeManager.getInstance().defineGenome(
                    genomeZipLocation, cytobandFileName, refFlatFileName,
                    fastaFileName, chrAliasFile, relativeSequenceLocation, genomeDisplayName,
                    genomeId, genomeFileName, monitor, seqLocationOverride);

            enableRemoveGenomes();

            igvCommandBar.addToUserDefinedGenomeItemList(genomeListItem);
            igvCommandBar.selectGenomeFromListWithNoImport(genomeListItem.getId());

            if (monitor != null) {
                monitor.fireProgressChange(100);
            }

        } catch (MaximumContigGenomeException e) {

            String genomePath = "";
            if (archiveFile != null) {
                genomePath = archiveFile.getAbsolutePath();
            }

            log.error("Failed to define genome: " + genomePath, e);

            JOptionPane.showMessageDialog(this, "Failed to define the current genome " +
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
        return igvCommandBar.getGenomeSelectedInDropdown();
    }

    /**
     * Gets the collection of genome display names currently in use.
     *
     * @return Set of display names.
     */
    public Collection<String> getGenomeDisplayNames() {
        return igvCommandBar.getGenomeDisplayNames();
    }

    public Collection<String> getGenomeIds() {
        return igvCommandBar.getGenomeIds();
    }

    public GenomeListItem doLoadGenome(ProgressMonitor monitor) {

        ProgressBar bar = null;
        GenomeListItem genomeListItem = null;
        boolean doImport = true;
        while (doImport) {

            doImport = false;
            File file = null;
            CursorToken token = WaitCursorManager.showWaitCursor();
            try {
                File importDirectory =
                        PreferenceManager.getInstance().getLastGenomeImportDirectory();
                if (importDirectory == null) {
                    PreferenceManager.getInstance().setLastGenomeImportDirectory(Globals.getUserDirectory());
                }

                FileFilter[] fileFilters = {new SnapshotUtilities.GenomeArchiveFileFilter()};

                genomeImportFileChooser = getFileChooser(importDirectory, null, fileFilters);
                genomeImportFileChooser.setDialogTitle("Load Genome");
                genomeImportFileChooser.addPropertyChangeListener(
                        new PropertyChangeListener() {

                            public void propertyChange(PropertyChangeEvent e) {

                                File oldFile = null;
                                String property = e.getPropertyName();
                                if (JFileChooser.SELECTED_FILE_CHANGED_PROPERTY.equals(property)) {
                                    oldFile = (File) e.getOldValue();
                                    genomeImportFileChooser.setPreviousFile(oldFile);
                                } else if (JFileChooser.FILE_FILTER_CHANGED_PROPERTY.equals(property)) {

                                    if (e.getOldValue() instanceof SnapshotUtilities.GenomeArchiveFileFilter &&
                                            e.getNewValue() instanceof SnapshotUtilities.GenomeArchiveFileFilter) {

                                        SnapshotUtilities.GenomeArchiveFileFilter newFilter =
                                                (SnapshotUtilities.GenomeArchiveFileFilter) e.getNewValue();

                                        File currentDirectory = genomeImportFileChooser.getCurrentDirectory();
                                        File previousFile = genomeImportFileChooser.getPreviousFile();
                                        if (previousFile != null) {

                                            File file = null;
                                            if (currentDirectory != null) {
                                                file = new File(currentDirectory,
                                                        previousFile.getName());
                                            } else {
                                                file = previousFile;
                                            }

                                            final File selectedFile = Utilities.changeFileExtension(
                                                    file, newFilter.getExtension());

                                            UIUtilities.invokeOnEventThread(new Runnable() {

                                                public void run() {
                                                    genomeImportFileChooser.setSelectedFile(
                                                            selectedFile);
                                                    genomeImportFileChooser.validate();
                                                }
                                            });
                                        }

                                    }
                                }
                            }
                        });

                // Display the dialog
                genomeImportFileChooser.showOpenDialog(this);
                file = genomeImportFileChooser.getSelectedFile();

                // If a file selection was made
                if (file != null) {
                    if (monitor != null) {
                        bar = ProgressBar.showProgressDialog(IGVMainFrame.getInstance(),
                                "Loading Genome...", monitor, false);
                    }

                    File directory = genomeImportFileChooser.getCurrentDirectory();
                    if (directory != null) {
                        PreferenceManager.getInstance().setLastGenomeImportDirectory(directory);
                    }

                    try {

                        if (monitor != null) {
                            monitor.fireProgressChange(50);
                        }

                        // Import the genome

                        if (log.isDebugEnabled()) {
                            log.debug("Call loadGenome");
                        }
                        genomeListItem = GenomeManager.getInstance().loadGenome(file.getAbsolutePath(), true, monitor);

                        igvCommandBar.addToUserDefinedGenomeItemList(genomeListItem);
                        igvCommandBar.selectGenomeFromListWithNoImport(genomeListItem.getId());


                        if (monitor != null) {
                            monitor.fireProgressChange(100);
                        }

                        if (bar != null) {
                            bar.close();
                        }

                    } catch (Exception e) {
                        log.fatal("Could not import genome!", e);
                    } finally {
                    }
                }
            } catch (Exception e) {

                String genomePath = "";
                if (file != null) {
                    genomePath = file.getAbsolutePath();
                }

                log.error("Failed to load genome: " + genomePath, e);
                int option =
                        JOptionPane.showConfirmDialog(this, "Failed to load the current genome " +
                                genomePath + "\n" + "Would you like to load another?",
                                "Load Genome Failure", JOptionPane.OK_CANCEL_OPTION);

                if (option == JOptionPane.OK_OPTION) {
                    doImport = true;
                }

            } finally {
                WaitCursorManager.removeWaitCursor(token);
            }

        }

        return genomeListItem;
    }

    private List<AbstractButton> createMenus() {

        List<AbstractButton> menus = new ArrayList<AbstractButton>();
        menus.add(createFileMenu());
        menus.add(createViewMenu());
        menus.add(createTracksMenu());

        extrasMenu = createExtrasMenu();
        extrasMenu.setVisible(false);
        menus.add(extrasMenu);

        menus.add(createHelpMenu());

        // Experimental -- remove for production release

        return menus;
    }

    public void enableExtrasMenu() {
        extrasMenu.setVisible(true);
    }

    /**
     * Load a collection of tracks in a background thread.
     *
     * @param locators
     */
    public void loadTracks(final Collection<ResourceLocator> locators) {
        loadTracks(locators, false);
    }


    /**
     * Load tracks corresponding to a collection of resource locations.
     * Note: Most of the code here is to adjust the scrollbars and split pane after loading
     * // TODO -- why is this in the batch frame (as opposed to TrackManager for example)?
     *
     * @param locators
     */
    public void loadTracks(final Collection<ResourceLocator> locators, boolean doInBackground) {

        ((ApplicationStatusBar) statusBar).setMessage("Loading ...");

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
                final JideSplitPane centerSplitPane = mainPanel.getCenterSplitPane();
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
                        mainPanel.doLayout();


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

    private JMenu createFileMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        menuItems.add(new JSeparator());

        // Load menu items
        menuAction = new LoadFilesMenuAction("Load from File...", KeyEvent.VK_L, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_FROM_URL, KeyEvent.VK_U, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new LoadFromServerAction("Load from Server...", KeyEvent.VK_S, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_FROM_DAS, KeyEvent.VK_D, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        // Session menu items
        menuAction = new NewSessionMenuAction("New Session...", KeyEvent.VK_N, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new OpenSessionMenuAction("Open Session...", KeyEvent.VK_O, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new SaveSessionMenuAction("Save Session...", KeyEvent.VK_V, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        menuAction =
                new MenuAction(UIConstants.IMPORT_GENOME_LIST_MENU_ITEM, null, KeyEvent.VK_D) {

                    @Override
                    public void actionPerformed(ActionEvent event) {

                        SwingWorker worker = new SwingWorker() {

                            public Object doInBackground() {

                                ProgressMonitor monitor = new ProgressMonitor();
                                doDefineGenome(monitor);
                                return null;
                            }
                        };
                        worker.execute();
                    }
                };

        menuAction.setToolTipText(IMPORT_GENOME_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        boolean hasImportedGenomes = true;
        try {
            hasImportedGenomes = !GenomeManager.getInstance().getUserDefinedGenomeArchiveList().isEmpty();

        } catch (IOException iOException) {
            // Ignore
        }
        removeImportedGenomeAction = new RemoveUserDefinedGenomeMenuAction(
                UIConstants.REMOVE_GENOME_LIST_MENU_ITEM, KeyEvent.VK_R);
        removeImportedGenomeAction.setEnabled(hasImportedGenomes);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(removeImportedGenomeAction));

        //menuAction = new ClearGenomeCacheAction("Clear Genome Cache...");
        //menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        // ***** Snapshots
        // Snapshot Application
        menuAction =
                new MenuAction("Save Image ...", null, KeyEvent.VK_A) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        doApplicationSnapshot(mainPanel);

                    }
                };

        menuAction.setToolTipText(SAVE_IMAGE_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        // Export Regions
        menuAction = new ExportRegionsMenuAction("Export Regions ...", KeyEvent.VK_E, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        // Import Regions
        menuAction = new ImportRegionsMenuAction("Import Regions ...", KeyEvent.VK_I, this);
        menuAction.setToolTipText(IMPORT_REGION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Import Regions
        menuAction = new ClearRegionsMenuAction("Clear Regions ...", this);
        menuAction.setToolTipText(IMPORT_REGION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        //dhmay adding 2010/11/16
        // Navigate Regions
        menuAction = new NavigateRegionsMenuAction("Navigate Regions ...", this);
        menuAction.setToolTipText(UIConstants.NAVIGATE_REGION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Separator
        /*menuItems.add(new JSeparator());
        menuAction =
                new MenuAction("Preprocess ...", null, KeyEvent.VK_P) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        (new PreprocessorDialog(IGVMainFrame.this, false)).setVisible(true);
                    }
                };

        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
        */

        // batch script
        menuItems.add(new JSeparator());

        menuAction = new RunScriptMenuAction("Run Batch Script...", KeyEvent.VK_X, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // igvtools
        menuItems.add(new JSeparator());
        menuAction = new SortTracksMenuAction("Compute coverage...", KeyEvent.VK_T, this) {
            @Override
            public void actionPerformed(ActionEvent e) {
                CoverageGui.launch(false, GenomeManager.getInstance().getGenomeId(), CoverageGui.Mode.COVERAGE);
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new SortTracksMenuAction("Convert to tdf...", KeyEvent.VK_T, this) {
            @Override
            public void actionPerformed(ActionEvent e) {
                CoverageGui.launch(false, GenomeManager.getInstance().getGenomeId(), CoverageGui.Mode.TILE);
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menuAction = new SortTracksMenuAction("Create index...", KeyEvent.VK_T, this) {
            @Override
            public void actionPerformed(ActionEvent e) {
                IndexGui.launch(false);
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new SortTracksMenuAction("Run igvtools...", KeyEvent.VK_T, this) {
            @Override
            public void actionPerformed(ActionEvent e) {
                IgvToolsGui.launch(false, GenomeManager.getInstance().getGenomeId());
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());      // Exit
        menuAction =
                new MenuAction("Exit", null, KeyEvent.VK_X) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        doExitApplication();
                    }
                };

        menuAction.setToolTipText(EXIT_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        // Empty the recent sessions list before we start to do
        // anything with it
        recentSessionList.clear();

        // Retrieve the stored session paths
        String recentSessions = PreferenceManager.getInstance().getRecentSessions();
        if (recentSessions != null) {
            String[] sessions = recentSessions.split(";");
            for (String sessionPath : sessions) {
                if (!recentSessionList.contains(sessionPath)) {
                    recentSessionList.add(sessionPath);
                }

            }
        }

        if (!recentSessionList.isEmpty()) {

            menuItems.add(new JSeparator());

            // Now add menu items
            for (final String session : recentSessionList) {
                OpenSessionMenuAction osMenuAction = new OpenSessionMenuAction(session, new File(session), this);
                menuItems.add(MenuAndToolbarUtils.createMenuItem(osMenuAction));
            }

        }

        MenuAction fileMenuAction = new MenuAction("File", null, KeyEvent.VK_F);
        return MenuAndToolbarUtils.createMenu(menuItems, fileMenuAction);
    }

    private JMenu createTracksMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        // Sort Context
        menuAction = new SortTracksMenuAction("Sort Tracks ...", KeyEvent.VK_S, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new GroupTracksMenuAction("Group Tracks  ... ", KeyEvent.VK_G, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Filter Tracks
        filterTracksAction = new FilterTracksMenuAction("Filter Tracks ...", KeyEvent.VK_F, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(filterTracksAction));

        menuItems.add(new JSeparator());

        // Reset Tracks
        menuAction = new FitDataToWindowMenuAction("Fit Data to Window", KeyEvent.VK_W, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        // Set track height
        menuAction = new SetTrackHeightMenuAction("Set Track Height...", KeyEvent.VK_H, this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        MenuAction dataMenuAction = new MenuAction("Tracks", null, KeyEvent.VK_K);
        return MenuAndToolbarUtils.createMenu(menuItems, dataMenuAction);
    }


    private JMenu createViewMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        // Preferences
        menuAction =
                new MenuAction("Preferences...", null, KeyEvent.VK_P) {

                    @Override
                    public void actionPerformed(ActionEvent e) {

                        UIUtilities.invokeOnEventThread(new Runnable() {

                            public void run() {
                                doViewPreferences();
                            }
                        });
                    }
                };
        menuAction.setToolTipText(PREFERENCE_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction =
                new MenuAction("Color Legends ...", null, KeyEvent.VK_H) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        (new LegendDialog(IGVMainFrame.this, false)).setVisible(true);
                    }
                };
        menuAction.setToolTipText(SHOW_HEATMAP_LEGEND_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        menuAction = new MenuAction("Show Name Panel", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                if (menuItem.isSelected()) {
                    mainPanel.expandNamePanel();
                } else {
                    mainPanel.collapseNamePanel();
                }
                doRefresh();
            }
        };
        boolean isShowing = mainPanel.isExpanded();
        JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
        menuItem.setSelected(isShowing);
        menuItem.setAction(menuAction);
        menuItems.add(menuItem);

        // Hide or Show the attribute panels
        boolean isShow = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY);
        doShowAttributeDisplay(isShow);  // <= WEIRD doing this here!

        menuAction = new MenuAction("Show Attribute Display", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                PreferenceManager.getInstance().setShowAttributeView(menuItem.getState());
                mainPanel.invalidate();
                doRefresh();
            }
        };
        menuAction.setToolTipText(SHOW_ATTRIBUTE_DISPLAY_TOOLTIP);
        menuItem = MenuAndToolbarUtils.createMenuItem(menuAction, isShow);
        menuItems.add(menuItem);


        menuAction =
                new MenuAction("Select Attributes to Show...", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        doSelectDisplayableAttribute();
                    }
                };
        menuAction.setToolTipText(SELECT_DISPLAYABLE_ATTRIBUTES_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        //
        menuAction =
                new MenuAction("Gene Lists...", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        (new GeneListManagerUI(IGVMainFrame.this)).setVisible(true);
                    }
                };
        menuAction.setToolTipText(SELECT_DISPLAYABLE_ATTRIBUTES_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        /*
        menuAction =
                new MenuAction("Show Region Bars", null, KeyEvent.VK_A) {

                    @Override
                    public void actionPerformed(ActionEvent e) {

                        JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                        PreferenceManager.getInstance().setShowRegionBars(menuItem.isSelected());
                        repaintDataPanels();
                    }
                };
                */
        //menuAction.setToolTipText(SHOW_ATTRIBUTE_DISPLAY_TOOLTIP);
        // menuItem = MenuAndToolbarUtils.createMenuItem(menuAction, PreferenceManager.getInstance().isShowRegionBars());
        // menuItems.add(menuItem);

        menuItems.add(new JSeparator());
        menuAction =
                new MenuAction("Refresh", null, KeyEvent.VK_R) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        doRefresh();
                    }
                };
        menuAction.setToolTipText(REFRESH_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new HistoryMenu("Go to"));


        // Add to IGVPanel menu
        MenuAction dataMenuAction = new MenuAction("View", null, KeyEvent.VK_V);
        return MenuAndToolbarUtils.createMenu(menuItems, dataMenuAction);
    }


    public void setGeneList(String listID) {
        setGeneList(listID, true);
    }

    public void setGeneList(String listID, boolean recordHistory) {
        GeneList gl = GeneListManager.getGeneList(listID);

        if (recordHistory) {
            session.getHistory().push("List: " + listID);
        }
        session.setCurrentGeneList(gl);
        resetFrames();

    }

    public void setDefaultFrame(String searchString) {
        FrameManager.setToDefaultFrame(searchString);
        resetFrames();
    }

    public void resetFrames() {
        mainPanel.headerPanelContainer.createHeaderPanels();
        for (TrackPanelScrollPane tp : trackManager.getTrackPanelScrollPanes()) {
            tp.getTrackPanel().createDataPanels();
        }

        igvCommandBar.setGeneListMode(FrameManager.isGeneListMode());
        mainPanel.revalidate();
        mainPanel.applicationHeaderPanel.revalidate();
        mainPanel.repaint();
    }


    private JMenu createExtrasMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();

        MenuAction menuAction = null;

        // Preferences reset
        menuAction = new ResetPreferencesAction("Reset Preferences", this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Linked sorting
        menuAction =
                new MenuAction("Use Linked Sorting", null, KeyEvent.VK_C) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        JCheckBoxMenuItem menuItem =
                                (JCheckBoxMenuItem) e.getSource();
                        boolean isSelected = menuItem.getState();
                        PreferenceManager.getInstance().put(PreferenceManager.ENABLE_LINKED_SORTING, String.valueOf(isSelected));
                    }
                };
        boolean linkedSorting = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_LINKED_SORTING);
        menuAction.setToolTipText("Enable linked sorting");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction, linkedSorting));


        menuItems.add(new JSeparator());

        // Load genome
        menuAction =
                new MenuAction(LOAD_GENOME_LIST_MENU_ITEM, null, KeyEvent.VK_I) {

                    @Override
                    public void actionPerformed(ActionEvent event) {

                        SwingWorker worker = new SwingWorker() {

                            public Object doInBackground() {
                                ProgressMonitor monitor = new ProgressMonitor();
                                doLoadGenome(monitor);
                                return null;
                            }
                        };
                        worker.execute();
                    }
                };

        menuAction.setToolTipText(LOAD_GENOME_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Set frame dimensions
        menuAction =
                new MenuAction("Set window dimensions", null, KeyEvent.VK_C) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        String value = JOptionPane.showInputDialog("Enter dimensions, e.g. 800x400");
                        String[] vals = value.split("x");
                        if (vals.length == 2) {
                            int w = Integer.parseInt(vals[0]);
                            int h = Integer.parseInt(vals[1]);
                            IGVMainFrame.getInstance().setSize(w, h);
                        }
                    }
                };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Save entire window
        menuAction =
                new MenuAction("Save Screenshot ...", null, KeyEvent.VK_A) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        doApplicationSnapshot(getContentPane());

                    }
                };

        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        //
        JMenu lfMenu = new JMenu("L&F");
        LookAndFeel lf = UIManager.getLookAndFeel();
        for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {

            final String lfName = info.getName();
            JMenuItem cb = new JMenuItem(lfName);
            //cb.setSelected(info.getClassName().equals(lf.getClass().getName());
            cb.addActionListener(new AbstractAction() {

                public void actionPerformed(ActionEvent actionEvent) {
                    for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {

                        if (lfName.equals(info.getName())) {
                            try {
                                UIManager.setLookAndFeel(info.getClassName());
                            } catch (ClassNotFoundException e) {
                                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                            } catch (InstantiationException e) {
                                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                            } catch (IllegalAccessException e) {
                                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                            } catch (UnsupportedLookAndFeelException e) {
                                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                            }
                            break;
                        }
                    }
                }
            });
            lfMenu.add(cb);
        }

        menuAction = new ExportTrackNamesMenuAction("Export track names...", this);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        MenuAction extrasMenuAction = new MenuAction("Extras");
        JMenu menu = MenuAndToolbarUtils.createMenu(menuItems, extrasMenuAction);

        menu.add(lfMenu);

        menu.setVisible(false);


        return menu;
    }

    private JMenu createHelpMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();

        MenuAction menuAction = null;

        menuAction =
                new MenuAction("Help ... ") {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        try {
                            BrowserLauncher.openURL(SERVER_BASE_URL + "igv/UserGuide");
                        } catch (IOException ex) {
                            log.error("Error opening browser", ex);
                        }

                    }
                };
        menuAction.setToolTipText(HELP_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menuAction =
                new MenuAction("Tutorial ... ") {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        try {
                            BrowserLauncher.openURL(SERVER_BASE_URL + "igv/QuickStart");
                        } catch (IOException ex) {
                            log.error("Error opening browser", ex);
                        }

                    }
                };
        menuAction.setToolTipText(TUTORIAL_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        if (Desktop.isDesktopSupported()) {
            final Desktop desktop = Desktop.getDesktop();
            if (desktop.isSupported(Desktop.Action.MAIL)) {

                menuAction =
                        new MenuAction("Contact Support") {

                            @Override
                            public void actionPerformed(ActionEvent e) {
                                try {
                                    URI uri = new URI("mailto:igv-help@broadinstitute.org");
                                    Desktop.getDesktop().mail(uri);
                                } catch (Exception ex) {
                                    log.error("Error opening email client", ex);
                                }

                            }
                        };
                menuAction.setToolTipText("Email support");
                menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
            }
        }

        menuAction =
                new MenuAction("About IGV ") {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        (new AboutDialog(IGVMainFrame.this, true)).setVisible(true);
                    }
                };
        menuAction.setToolTipText(ABOUT_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        MenuAction toolMenuAction = new MenuAction("Help");
        return MenuAndToolbarUtils.createMenu(menuItems, toolMenuAction);
    }

    public void enableRemoveGenomes() {
        if (removeImportedGenomeAction != null) {
            removeImportedGenomeAction.setEnabled(true);
        }

    }


    /**
     * Select a genome
     */
    final public void doChooseGenome(GenomeDescriptor genomeType) {

        CursorToken token = null;
        try {

            token = WaitCursorManager.showWaitCursor();

            if (genomeType != null) {

                final String genomeId = genomeType.getId();
                String currentGenomeId = GenomeManager.getInstance().getGenomeId();
                if (currentGenomeId != null && genomeId.equalsIgnoreCase(currentGenomeId)) {
                    // Nothing to do if genome already loaded
                    return;
                }

                setGenomeId(genomeId);
                PreferenceManager.getInstance().setDefaultGenome(genomeId);
                IGVMainFrame.getInstance().getTrackManager().reloadSAMTracks();
            }

        } finally {
            WaitCursorManager.removeWaitCursor(token);
        }

    }

    /**
     * Open the user preferences dialog
     */
    final public void doViewPreferences() {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                boolean originalSingleTrackValue =
                        PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY);

                PreferencesEditor dialog = new PreferencesEditor(IGVMainFrame.this, true);
                dialog.setVisible(true);


                if (dialog.isCanceled()) {
                    resetStatusMessage();
                    return;

                }


                try {

                    //Should data and feature panels be combined ?
                    boolean singlePanel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY);
                    if (originalSingleTrackValue != singlePanel) {
                        JOptionPane.showMessageDialog(IGVMainFrame.this, "Panel option change will take affect after restart.");
                    }


                } finally {

                    // Update the state of the current tracks for drawing purposes
                    updateTrackState();
                    resetStatusMessage();

                }


            }
        });
    }

    final public void doExitApplication() {

        try {

            ((ApplicationStatusBar) statusBar).setMessage("Exiting...");

            // Store recent sessions
            if (!recentSessionList.isEmpty()) {

                int size = recentSessionList.size();
                if (size > UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST) {
                    size = UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST;
                }

                String recentSessions = "";
                for (int i = 0; i <
                        size; i++) {
                    recentSessions += recentSessionList.get(i);

                    if (i < (size - 1)) {
                        recentSessions += ";";
                    }

                }
                PreferenceManager.getInstance().remove(PreferenceManager.RECENT_SESSION_KEY);
                PreferenceManager.getInstance().setRecentSessions(recentSessions);
            }

            PreferenceManager.getInstance().setApplicationFrameBounds(getBounds());

            setVisible(false);
        } finally {
            System.exit(0);
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

        rootPane.revalidate();
        repaint();
        //getContentPane().repaint();
    }

    final public void refreshCommandBar() {
        igvCommandBar.updateCurrentCoordinates();
    }


// TODO -- move all of this attribute stuf out of IGVMainFrame,  perhaps to

    // some Attribute helper class.

    final public void doSelectDisplayableAttribute() {

        List<String> allAttributes = AttributeManager.getInstance().getAttributeKeys();
        Set<String> hiddenAttributes = AttributeManager.getInstance().getHiddenAttributes();
        final CheckListDialog dlg = new CheckListDialog(this, allAttributes, hiddenAttributes, false);
        dlg.setVisible(true);

        if (!dlg.isCanceled()) {
            AttributeManager.getInstance().setHiddenAttributes(dlg.getNonSelections());
            doRefresh();
        }
    }


    final public void doApplicationSnapshot(Component target) {
        ((ApplicationStatusBar) statusBar).setMessage("Creating snapshot...");
        File defaultFile = new File("igv_snapshot.png");
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
            ((ApplicationStatusBar) statusBar).setMessage("Exporting image: " + defaultFile.getAbsolutePath());
            File file = selectSnapshotFile(defaultFile);
            if (file == null) {
                return;
            }
            createSnapshotNonInteractive(target, file);
        } catch (Exception e) {
            log.error("Error creating exporting image ", e);
            MessageUtils.showMessage(("Error creating the image file: " + defaultFile + "<br> "
                    + e.getMessage()));
        }
        finally {
            WaitCursorManager.removeWaitCursor(token);
            resetStatusMessage();
        }

    }


    public void createSnapshotNonInteractive(File file) {
        createSnapshotNonInteractive(mainPanel, file);
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

            boolean doubleBuffered = RepaintManager.currentManager(getContentPane()).isDoubleBufferingEnabled();
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
        snapshotFileChooser.showSaveDialog(this);

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

    public void setGenomeId(String id) {

        if (log.isDebugEnabled()) {
            log.debug("Setting current genome id");
        }

        String currentGenomeId = GenomeManager.getInstance().getGenomeId();
        if (currentGenomeId != null && id.equalsIgnoreCase(currentGenomeId)) {
            // Nothing to do if genome already loaded
            return;
        }

        String gid = GenomeManager.getInstance().setGenomeId(id);

        FeatureDB.clearFeatures();

        IGVMainFrame.getInstance().getTrackManager().loadGeneTrack(gid);


        for (Chromosome chr : GenomeManager.getInstance().getCurrentGenome().getChromosomes()) {
            for (Cytoband cyto : chr.getCytobands()) {
                FeatureDB.addFeature(cyto.getLongName(), cyto);
            }
        }


        if (igvCommandBar != null) {
            igvCommandBar.updateChromosomeDropdown();
        }

        PreferenceManager.getInstance().setDefaultGenome(gid);
    }

    private void createZoomCursors() throws HeadlessException, IndexOutOfBoundsException {
        if (zoomInCursor == null || zoomOutCursor == null) {
            final Image zoomInImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_IN).getImage();
            final Image zoomOutImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_OUT).getImage();
            final Point hotspot = new Point(10, 10);
            zoomInCursor = getToolkit().createCustomCursor(zoomInImage, hotspot, "Zoom in");
            zoomOutCursor = getToolkit().createCustomCursor(zoomOutImage, hotspot, "Zoom out");

        }

    }

    private void createHandCursor() throws HeadlessException, IndexOutOfBoundsException {
        if (handCursor == null) {
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
        }

        if (fistCursor == null) {
            BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

            // Make backgroun transparent
            Graphics2D g = handImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(
                    AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
            g.fill(rect);

            // Draw hand image in middle
            g =
                    handImage.createGraphics();
            g.drawImage(IconFactory.getInstance().getIcon(
                    IconFactory.IconID.FIST).getImage(), 0, 0, null);
            fistCursor =
                    getToolkit().createCustomCursor(
                            handImage, new Point(8, 6), "Move");
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
                    getToolkit().createCustomCursor(
                            dragNDropImage, new Point(0, 0), "Drag and Drop");
        }

    }

    public void createNewSession(String sessionName) {

        LRUCache.clearCaches();

        AttributeManager.getInstance().clearAllAttributes();

        setTitle(UIConstants.APPLICATION_NAME);

        if (filterTracksAction != null) {
            filterTracksAction.resetTrackFilter();
        }

        AttributeManager.getInstance().clearAllAttributes();
        session = new Session(sessionName);

        mainPanel.resetPanels();


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

        ((ApplicationStatusBar) statusBar).setMessage(message);
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
            MessageUtils.showAndLogErrorMessage(IGVMainFrame.theInstance, message, log, e);
        }

    }

    public void updateTrackState() {

        doRefresh();
    }


    public void updateTrackFilter() {
        if (filterTracksAction != null) {
            // TODO -- HOLDING THIS STATE ON THE ACTION OBJECT SEEMS BAD
            filterTracksAction.updateTrackFilter();
            // Update the state of the current tracks
            updateTrackState();
        }
    }

    public void setFilterMatchAll(boolean value) {
        if (filterTracksAction != null) {
            filterTracksAction.setFilterMatchAll(value);
        }

    }

    public boolean isFilterMatchAll() {
        if (filterTracksAction != null) {
            return filterTracksAction.isFilterMatchAll();
        }

        return false;
    }

    public void setFilterShowAllTracks(boolean value) {
        if (filterTracksAction != null) {
            filterTracksAction.setFilterShowAllTracks(value);
        }

    }

    public boolean isFilterShowAllTracks() {
        if (filterTracksAction != null) {
            return filterTracksAction.getShowAllTracksFilterCheckBox().isSelected();
        }

        return false;
    }

    /**
     * Add a new data panel set
     */
    public TrackPanelScrollPane addDataPanel(String name) {

        return mainPanel.addDataPanel(name);
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
        return mainPanel.getTrackPanels();
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
                if (!recentSessionList.contains(sessionFilePath)) {
                    recentSessionList.addFirst(sessionFilePath);
                }

            } catch (Exception e) {
                String message = "Failed to load session! : " + sessionFile.getAbsolutePath();
                MessageUtils.showAndLogErrorMessage(IGVMainFrame.this, message, log, e);
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
            MessageUtils.showAndLogErrorMessage(IGVMainFrame.this, message, log);
        }

    }

    /**
     * TODO -- this is nearly an exact copy of the doRestoreSession(File sessionFile)
     * method.  Refactor to combine these using streams.
     *
     * @param sessionURL
     */
    final public void doRestoreSession(final URL sessionURL,
                                       final String locus) {

        if (log.isDebugEnabled()) {
            log.debug("Enter doRestoreSession: " + sessionURL + " " + locus);
        }

        if (sessionURL != null) {
            InputStream inputStream = null;
            try {
                inputStream = new BufferedInputStream(sessionURL.openStream());
                doRestoreSession(inputStream, URLDecoder.decode(sessionURL.getFile(), "UTF-8"), locus, false);
            } catch (Exception e) {
                String message = "Failed to load session! : " + sessionURL;
                MessageUtils.showAndLogErrorMessage(IGVMainFrame.this, message, log, e);
            } finally {

                if (inputStream != null) {
                    try {
                        inputStream.close();
                    } catch (IOException iOException) {
                        log.error("Error closing session stream", iOException);
                    }

                }
            }


        } else {
            String message = "Session file does not exist! : ";
            try {
                message += URLDecoder.decode(sessionURL.getFile(), "UTF-8");
            } catch (UnsupportedEncodingException ex) {
                message += sessionURL.getFile();
            }

            MessageUtils.showAndLogErrorMessage(IGVMainFrame.this, message, log);
        }

        if (log.isDebugEnabled()) {
            log.debug("Exit doRestoreSession");
        }

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

            (new SessionReader()).loadSession(inputStream, session, sessionPath);
            String searchText = locus == null ? session.getLocusString() : locus;

            // NOTE: Nothing to do if chr == all
            if (!FrameManager.isGeneListMode() && searchText != null && !searchText.equals(Globals.CHR_ALL) && searchText.trim().length() > 0) {
                igvCommandBar.searchByLocus(searchText);


            }

            setTitle(UIConstants.APPLICATION_NAME + " - Session: " + sessionPath);
            LRUCache.clearCaches();
            doRefresh();

            //If there's a RegionNavigatorDialog, kill it.
            //this could be done through the Observer that RND uses, I suppose.  Not sure that's cleaner
            RegionNavigatorDialog.destroyActiveInstance();
        } catch (Exception e) {
            String message = "Failed to load session! : " + sessionPath;
            MessageUtils.showAndLogErrorMessage(IGVMainFrame.this, message, log, e);
        } finally {

            resetStatusMessage();
        }

    }

    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        ((ApplicationStatusBar) statusBar).setMessage("" +
                IGVMainFrame.getInstance().getTrackManager().getVisibleTrackCount() + " tracks loaded");

    }


    public void rebuildGenomeDropdownList(Set excludedArchivesUrls) {
        igvCommandBar.rebuildGenomeItemList(excludedArchivesUrls);
    }

    public void showLoadedTrackCount() {
        ((ApplicationStatusBar) statusBar).setMessage("" +
                IGVMainFrame.getInstance().getTrackManager().getVisibleTrackCount() +
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

        igvCommandBar.searchByLocus(locus);
    }


    public TrackManager getTrackManager() {
        return trackManager;
    }

    public void tweakPanelDivider() {
        mainPanel.tweakPanelDivider();
    }

    public void removeDataPanel(String name) {
        mainPanel.removeDataPanel(name);
    }

    public void layoutMainPanel() {
        mainPanel.doLayout();
    }

    public Component getMainPanel() {
        return mainPanel;
    }

    public void setExportingSnapshot(boolean exportingSnapshot) {
        isExportingSnapshot = exportingSnapshot;
        if (isExportingSnapshot) {
            RepaintManager.currentManager(getContentPane()).setDoubleBufferingEnabled(false);
        } else {
            RepaintManager.currentManager(getContentPane()).setDoubleBufferingEnabled(true);
        }
    }


    /**
     * Startup the IGV batch window,  then execute batch file if supplied.
     *
     * @param args
     */
    public void startUp(final String[] args) {

        if (log.isDebugEnabled()) {
            log.debug("startUp");
        }

        Main.IGVArgs igvArgs = new Main.IGVArgs(args);
        SwingWorker worker = new StartupWorker(igvArgs);
        worker.execute();


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
            final ProgressBar progressBar =
                    ProgressBar.showProgressDialog(IGVMainFrame.this, "Initializing Genome...", monitor, false);
            monitor.fireProgressChange(10);

            // Load the last genome and chromosome
            String genomeId = PreferenceManager.getInstance().getDefaultGenome();
            if (genomeId == null) {
                genomeId = GenomeManager.getInstance().getTopGenomeListItem().getId();
            }
            setGenomeId(genomeId);
            monitor.fireProgressChange(50);

            genomeId = GenomeManager.getInstance().getGenomeId(); // <= might have changed
            try {
                igvCommandBar.initializeGenomeList(monitor);
                igvCommandBar.selectGenomeFromListWithNoImport(genomeId);
            } catch (FileNotFoundException ex) {
                JOptionPane.showMessageDialog(IGVMainFrame.this, "Error initializing genome list: " + ex.getMessage());
                log.error("Error initializing genome list: ", ex);
            } catch (NoRouteToHostException ex) {
                JOptionPane.showMessageDialog(IGVMainFrame.this, "Network error initializing genome list: " + ex.getMessage());
                log.error("Network error initializing genome list: ", ex);
            }

            // Done
            closeWindow(progressBar);

            // Optional arguments
            if (igvArgs.getPropertyFile() != null) {

            }
            if (igvArgs.getDataServerURL() != null) {
                PreferenceManager.getInstance().overrideDataServerURL(igvArgs.getDataServerURL());
            }
            if (igvArgs.getGenomeServerURL() != null) {
                PreferenceManager.getInstance().overrideGenomeServerURL(igvArgs.getGenomeServerURL());
            }

            // Start up a port listener.  Port # can be overriden with "-p" command line switch
            boolean portEnabled = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.PORT_ENABLED);
            String portString = igvArgs.getPort();
            if (portEnabled || portString != null) {
                // Command listner thread
                int port = PreferenceManager.getInstance().getAsInt(PreferenceManager.PORT_NUMBER);
                if (portString != null) {
                    port = Integer.parseInt(portString);
                }
                CommandListener.start(port);
            }

            //If there is an argument assume it is a session file or url
            if (igvArgs.getSessionFile() != null || igvArgs.getDataFileString() != null) {

                if (log.isDebugEnabled()) {
                    log.debug("Loadding session data");
                }

                final IndefiniteProgressMonitor indefMonitor = new IndefiniteProgressMonitor(60);
                final ProgressBar bar2 = ProgressBar.showProgressDialog(IGVMainFrame.this, "Loading session data", indefMonitor, false);

                int idx = 0;


                indefMonitor.start();
                try {

                    if (log.isDebugEnabled()) {
                        log.debug("Calling restore session");
                    }


                    if (igvArgs.getGenomeId() != null) {
                        selectGenomeFromList(igvArgs.getGenomeId());
                    }


                    if (igvArgs.getSessionFile() != null) {
                        if (IGVHttpUtils.isURL(igvArgs.getSessionFile())) {
                            URL url = new URL(igvArgs.getSessionFile());
                            doRestoreSession(url, igvArgs.getLocusString());
                        } else {
                            File sf = new File(igvArgs.getSessionFile());
                            if (sf.exists()) {
                                doRestoreSession(sf, igvArgs.getLocusString());
                            }
                        }
                        doRefresh();
                    } else if (igvArgs.getDataFileString() != null) {
                        // Not an xml file, assume its a list of data files
                        String[] tokens = igvArgs.getDataFileString().split(",");
                        List<ResourceLocator> locators = new ArrayList();
                        for (String p : tokens) {
                            locators.add(new ResourceLocator(p));
                        }
                        getTrackManager().loadResources(locators);
                        doRefresh();
                    }


                } catch (Exception ex) {
                    String tmp = igvArgs.getSessionFile() != null ? igvArgs.getSessionFile() : igvArgs.getDataFileString();
                    JOptionPane.showMessageDialog(IGVMainFrame.this, "<html>Error loading session: " + tmp + "<br>" + ex.toString());
                    log.error("Error loading session: " + tmp, ex);
                }


                indefMonitor.stop();
                closeWindow(bar2);
            }

            if (igvArgs.getLocusString() != null) {
                goToLocus(igvArgs.getLocusString());
            }


            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    setVisible(true);
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

    // THis is her for backward compatibility
    public static void main(String [] args) {
        Main.main(args);
    }
}
