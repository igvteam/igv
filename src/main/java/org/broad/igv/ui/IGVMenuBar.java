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

package org.broad.igv.ui;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.aws.S3LoadDialog;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.charts.ScatterPlotUtils;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.commandbar.RemoveGenomesDialog;
import org.broad.igv.google.GoogleUtils;
import org.broad.igv.google.OAuthProvider;
import org.broad.igv.google.OAuthUtils;
import org.broad.igv.lists.GeneListManagerUI;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.tools.IgvToolsGui;
import org.broad.igv.tools.motiffinder.MotifFinderPlugin;
import org.broad.igv.track.CombinedDataSourceDialog;
import org.broad.igv.ui.action.*;
import org.broad.igv.ui.commandbar.GenomeComboBox;
import org.broad.igv.ui.legend.LegendDialog;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.ReorderPanelsDialog;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.encode.EncodeFileBrowser;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import javax.swing.plaf.basic.BasicBorders;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static org.broad.igv.prefs.Constants.*;
import static org.broad.igv.ui.UIConstants.*;

/**
 * Main menu bar at top of window. File / genomes / view / etc.
 * Singleton
 *
 * @author jrobinso
 * @date Apr 4, 2011
 */
public class IGVMenuBar extends JMenuBar implements IGVEventObserver {

    private static Logger log = Logger.getLogger(IGVMenuBar.class);
    public static final String GENOMESPACE_REG_TOOLTIP = "Register for GenomeSpace";
    public static final String GENOMESPACE_REG_PAGE = "http://www.genomespace.org/register";

    private JMenu extrasMenu;
    private FilterTracksMenuAction filterTracksAction;
    private JMenu viewMenu;
    IGV igv;

    private JMenu toolsMenu;

    /**
     * We store this as a field because we alter it if
     * we can't access genome server list
     */
    private JMenuItem loadFromServerMenuItem;

    private static final String LOAD_GENOME_SERVER_TOOLTIP = "Select genomes available on the server to appear in menu.";
    private static final String CANNOT_LOAD_GENOME_SERVER_TOOLTIP = "Could not reach genome server";

    private static IGVMenuBar instance;
    private JMenu googleMenu;
    private JMenu AWSMenu;
    private JMenuItem encodeMenuItem;
    private JMenuItem reloadSessionItem;

    public void notifyGenomeServerReachable(boolean reachable) {
        if (loadFromServerMenuItem != null) {
            UIUtilities.invokeOnEventThread(() -> {
                loadFromServerMenuItem.setEnabled(reachable);
                String tooltip = reachable ? LOAD_GENOME_SERVER_TOOLTIP : CANNOT_LOAD_GENOME_SERVER_TOOLTIP;
                loadFromServerMenuItem.setToolTipText(tooltip);
            });
        }
    }

    public void showAboutDialog() {
        (new AboutDialog(IGV.getMainFrame(), true)).setVisible(true);
    }

    static IGVMenuBar createInstance(IGV igv) {
        if (instance != null) {
            if (igv == instance.igv) {
                return instance;
            }
            throw new IllegalStateException("Cannot create another IGVMenuBar, use getInstance");
        }
        UIUtilities.invokeAndWaitOnEventThread(() -> instance = new IGVMenuBar(igv));
        return instance;
    }

    public static IGVMenuBar getInstance() {
        return instance;
    }

    private IGVMenuBar(IGV igv) {
        this.igv = igv;
        setBorder(new BasicBorders.MenuBarBorder(Color.GRAY, Color.GRAY));
        setBorderPainted(true);

        for (AbstractButton menu : createMenus()) {
            add(menu);
        }

        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, this);


        //This is for Macs, so showing the about dialog
        //from the command bar does what we want.
        if (Globals.IS_MAC) {
            DesktopIntegration.setAboutHandler(this);
            DesktopIntegration.setQuitHandler();
        }
    }

    private List<AbstractButton> createMenus() {

        List<AbstractButton> menus = new ArrayList<AbstractButton>();

        menus.add(createFileMenu());
        menus.add(createGenomesMenu());
        menus.add(createViewMenu());
        menus.add(createTracksMenu());
        menus.add(createRegionsMenu());

        refreshToolsMenu();
        menus.add(toolsMenu);

        extrasMenu = createExtrasMenu();
        //extrasMenu.setVisible(false);
        menus.add(extrasMenu);

        try {
            googleMenu = createGoogleMenu();
            if (googleMenu != null) {
                googleMenu.setVisible(PreferencesManager.getPreferences().getAsBoolean(ENABLE_GOOGLE_MENU));
                menus.add(googleMenu);
            }
        } catch (IOException e) {
            log.error("Error creating google menu: " + e.getMessage());
        }

        try {
            AWSMenu = createAWSMenu();
            AWSMenu.setVisible(AmazonUtils.isAWSProviderPresent());
            menus.add(AWSMenu);

        } catch (IOException e) {
            log.error("Error creating the Amazon AWS menu: " + e.getMessage());
            AWSMenu.setVisible(false);
        }

        menus.add(createHelpMenu());

        // Experimental -- remove for production release

        return menus;
    }

    public void updateAWSMenu() {
        AWSMenu.setVisible(AmazonUtils.isAWSProviderPresent());
    }

    /**
     * Generate the "tools" menu.
     * This is imperative, it is written to field {@code toolsMenu}.
     * Reason being, when we add (TODO remove)
     * a new tool, we need to refresh just this menu
     */
    void refreshToolsMenu() {
        List<JComponent> menuItems = new ArrayList<JComponent>(10);

        // batch script
        MenuAction menuAction = new RunScriptMenuAction("Run Batch Script...", KeyEvent.VK_X, igv);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // igvtools
        menuAction = new SortTracksMenuAction("Run igvtools...", KeyEvent.VK_T, igv) {
            @Override
            public void actionPerformed(ActionEvent e) {
                IgvToolsGui.launch(false, igv.getGenomeManager().getGenomeId());
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Motif finder
        menuItems.add(MotifFinderPlugin.getMenuItem());

        // BLAT
        menuItems.add(BlatClient.getMenuItem());

        // Combine data tracks
        JMenuItem combineDataItem = new JMenuItem("Combine Data Tracks");
        combineDataItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                CombinedDataSourceDialog dialog = new CombinedDataSourceDialog(IGV.getMainFrame());
                dialog.setVisible(true);
            }
        });
        menuItems.add(combineDataItem);


        MenuAction toolsMenuAction = new MenuAction("Tools", null);
        if (toolsMenu == null) {
            toolsMenu = MenuAndToolbarUtils.createMenu(menuItems, toolsMenuAction);
            toolsMenu.setName("Tools");
        } else {
            toolsMenu.removeAll();
            for (JComponent item : menuItems) {
                toolsMenu.add(item);
            }
        }

    }

    public void enableExtrasMenu() {
        extrasMenu.setVisible(true);
    }


    JMenu createFileMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        menuItems.add(new JSeparator());

        // Load menu items
        menuAction = new LoadFilesMenuAction("Load from File...", KeyEvent.VK_L, igv);
        menuAction.setToolTipText(UIConstants.LOAD_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_FROM_URL, KeyEvent.VK_U, igv);
        menuAction.setToolTipText(UIConstants.LOAD_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new LoadFromServerAction("Load from Server...", KeyEvent.VK_S, igv);
        menuAction.setToolTipText(UIConstants.LOAD_SERVER_DATA_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

//        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_FROM_DAS, KeyEvent.VK_D, igv);
//        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        if (PreferencesManager.getPreferences().getAsBoolean(DB_ENABLED)) {
            menuAction = new LoadFromDatabaseAction("Load from Database...", 0, igv);
            menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
        }


        encodeMenuItem = MenuAndToolbarUtils.createMenuItem(new BrowseEncodeAction("Load from ENCODE (2012)...", KeyEvent.VK_E, igv));
        menuItems.add(encodeMenuItem);
        String genomeId = IGV.getInstance().getGenomeManager().getGenomeId();
        encodeMenuItem.setVisible(EncodeFileBrowser.genomeSupported(genomeId));


        menuItems.add(new JSeparator());
        menuAction = new ReloadTracksMenuAction("Reload Tracks", -1, igv);
        menuAction.setToolTipText(RELOAD_SESSION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        // Session menu items
        menuAction = new NewSessionMenuAction("New Session...", KeyEvent.VK_N, igv);
        menuAction.setToolTipText(UIConstants.NEW_SESSION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new OpenSessionMenuAction("Open Session...", KeyEvent.VK_O, igv);
        menuAction.setToolTipText(OPEN_SESSION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new SaveSessionMenuAction("Save Session...", KeyEvent.VK_V, igv);
        menuAction.setToolTipText(UIConstants.SAVE_SESSION_TOOLTIP);
        JMenuItem saveSessionItem = MenuAndToolbarUtils.createMenuItem(menuAction);
        menuItems.add(saveSessionItem);

        menuAction = new ReloadSessionMenuAction("Reload Session", -1, igv);
        menuAction.setToolTipText(RELOAD_SESSION_TOOLTIP);
        reloadSessionItem = MenuAndToolbarUtils.createMenuItem(menuAction);
        reloadSessionItem.setEnabled(false);
        menuItems.add(reloadSessionItem);

        menuItems.add(new JSeparator());

        // ***** Snapshots
        // Snapshot Application
        menuAction =
                new MenuAction("Save Image ...", null, KeyEvent.VK_A) {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        igv.saveImage(igv.getMainPanel());

                    }
                };

        menuAction.setToolTipText(SAVE_IMAGE_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // TODO -- change "Exit" to "Close" for BioClipse
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
        igv.getRecentSessionList().clear();

        // Retrieve the stored session paths
        String recentSessions = PreferencesManager.getPreferences().getRecentSessions();
        if (recentSessions != null) {
            String[] sessions = recentSessions.split(";");
            for (String sessionPath : sessions) {
                if (!sessionPath.equals("null") &&
                        !igv.getRecentSessionList().contains(sessionPath)) {
                    igv.getRecentSessionList().add(sessionPath);
                }

            }
        }

        if (!IGV.getInstance().getRecentSessionList().isEmpty()) {
            menuItems.add(new JSeparator());
            // Now add menu items
            for (final String session : IGV.getInstance().getRecentSessionList()) {
                OpenSessionMenuAction osMenuAction = new OpenSessionMenuAction(session, session, IGV.getInstance());
                menuItems.add(MenuAndToolbarUtils.createMenuItem(osMenuAction));
            }

        }

        MenuAction fileMenuAction = new MenuAction("File", null, KeyEvent.VK_F);
        JMenu fileMenu = MenuAndToolbarUtils.createMenu(menuItems, fileMenuAction);

        return fileMenu;
    }

    private JMenu createGenomesMenu() {
        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        // Load genome
        menuAction =
                new MenuAction("Load Genome from File...", null, KeyEvent.VK_I) {
                    @Override
                    public void actionPerformed(ActionEvent event) {
                        try {
                            File importDirectory = PreferencesManager.getPreferences().getLastGenomeImportDirectory();
                            if (importDirectory == null) {
                                PreferencesManager.getPreferences().setLastGenomeImportDirectory(DirectoryManager.getUserDirectory());
                            }
                            // Display the dialog
                            File file = FileDialogUtils.chooseFile("Load Genome", importDirectory, FileDialog.LOAD);

                            // If a file selection was made
                            if (file != null) {
                                GenomeManager.getInstance().loadGenome(file.getAbsolutePath(), null);
                            }
                        } catch (Exception e) {
                            MessageUtils.showErrorMessage(e.getMessage(), e);
                        }
                    }
                };

        menuAction.setToolTipText("Load a FASTA, .json, or .genome file...");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Load genome from URL
        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_GENOME_FROM_URL, 0, igv);
        menuAction.setToolTipText("Load a FASTA, .json, or .genome file...");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());
        // Download genome from server 
        menuAction = new MenuAction("Download Genomes...", null) {
            @Override
            public void actionPerformed(ActionEvent event) {
                GenomeComboBox.loadGenomeFromServer();
            }
        };
        menuAction.setToolTipText(LOAD_GENOME_SERVER_TOOLTIP);
        loadFromServerMenuItem = MenuAndToolbarUtils.createMenuItem(menuAction);
        menuItems.add(loadFromServerMenuItem);

        menuItems.add(new JSeparator());

        // Add genome to combo box from server
        menuAction = new MenuAction("Remove Genomes...", null) {
            @Override
            public void actionPerformed(ActionEvent event) {
                RemoveGenomesDialog dialog2 = new RemoveGenomesDialog(IGV.getMainFrame());
                dialog2.setVisible(true);
            }
        };
        menuAction.setToolTipText("Remove genomes which appear in the dropdown list");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menuItems.add(new JSeparator());

        menuAction =
                new MenuAction("Create .genome File...", null, KeyEvent.VK_D) {
                    @Override
                    public void actionPerformed(ActionEvent event) {
                        javax.swing.ProgressMonitor monitor = new javax.swing.ProgressMonitor(IGV.getInstance().getMainPanel(),
                                "Creating genome", null, 0, 100);
                        igv.defineGenome(monitor);
                    }
                };

        menuAction.setToolTipText(UIConstants.IMPORT_GENOME_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        MenuAction genomeMenuAction = new MenuAction("Genomes", null);
        return MenuAndToolbarUtils.createMenu(menuItems, genomeMenuAction);
    }


    private JMenu createTracksMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        // Sort Context
        menuAction = new SortTracksMenuAction("Sort Tracks...", KeyEvent.VK_S, IGV.getInstance());
        menuAction.setToolTipText(SORT_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new GroupTracksMenuAction("Group Tracks... ", KeyEvent.VK_G, IGV.getInstance());
        menuAction.setToolTipText(UIConstants.GROUP_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Filter Tracks
        filterTracksAction = new FilterTracksMenuAction("Filter Tracks...", KeyEvent.VK_F, IGV.getInstance());
        filterTracksAction.setToolTipText(UIConstants.FILTER_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(filterTracksAction));

        // Rename tracks
        menuAction = new RenameTracksMenuAction("Rename Tracks... ", KeyEvent.VK_R, IGV.getInstance());
        menuAction.setToolTipText(UIConstants.RENAME_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Overlay tracks
        menuAction = new OverlayTracksMenuAction("Overlay Data Tracks... ", KeyEvent.VK_O, IGV.getInstance());
        menuAction.setToolTipText(UIConstants.OVERLAY_TRACKS_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menuItems.add(new JSeparator());

        // Reset Tracks
        menuAction = new FitDataToWindowMenuAction("Fit Data to Window", KeyEvent.VK_W, IGV.getInstance());
        menuAction.setToolTipText(UIConstants.FIT_DATA_TO_WINDOW_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        // Set track height
        menuAction = new SetTrackHeightMenuAction("Set Track Height...", KeyEvent.VK_H, IGV.getInstance());
        menuAction.setToolTipText(UIConstants.SET_DEFAULT_TRACK_HEIGHT_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        MenuAction dataMenuAction = new MenuAction("Tracks", null, KeyEvent.VK_K);

        //menuItems.add(exportData);

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
                                IGV.getInstance().doViewPreferences();
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
                        (new LegendDialog(IGV.getMainFrame())).setVisible(true);
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
                    IGV.getInstance().getMainPanel().expandNamePanel();
                } else {
                    IGV.getInstance().getMainPanel().collapseNamePanel();
                }
                IGV.getInstance().repaint();
            }
        };
        boolean isShowing = IGV.getInstance().getMainPanel().isExpanded();
        JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
        menuItem.setSelected(isShowing);
        menuItem.setAction(menuAction);
        menuItems.add(menuItem);

        JMenuItem panelWidthmenuItem = new JMenuItem();
        menuAction = new MenuAction("Set Name Panel Width...", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {
                MainPanel mainPanel = IGV.getInstance().getMainPanel();
                String currentValue = String.valueOf(mainPanel.getNamePanelWidth());
                String newValue = MessageUtils.showInputDialog("Enter track name panel width: ", currentValue);
                if (newValue != null) {
                    try {
                        Integer w = Integer.parseInt(newValue);
                        if (w <= 0 || w == 1000) throw new NumberFormatException();
                        PreferencesManager.getPreferences().put(NAME_PANEL_WIDTH, newValue);
                        mainPanel.setNamePanelWidth(w);
                    } catch (NumberFormatException ex) {
                        MessageUtils.showErrorMessage("Error: value must be a positive integer < 1000.", ex);
                    }
                }
            }
        };
        panelWidthmenuItem.setAction(menuAction);
        menuItems.add(panelWidthmenuItem);

        // Hide or Show the attribute panels
        boolean isShow = PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY);
        IGV.getInstance().doShowAttributeDisplay(isShow);  // <= WEIRD doing IGV.getInstance() here!

        menuAction = new MenuAction("Show Attribute Display", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                PreferencesManager.getPreferences().setShowAttributeView(menuItem.getState());
                IGV.getInstance().getMainPanel().invalidate();
                IGV.getInstance().repaint();
            }
        };
        menuItem = MenuAndToolbarUtils.createMenuItem(menuAction, isShow);
        menuItems.add(menuItem);


        menuAction =
                new MenuAction("Select Attributes to Show...", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        IGV.getInstance().doSelectDisplayableAttribute();
                    }
                };
        menuAction.setToolTipText(SELECT_DISPLAYABLE_ATTRIBUTES_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new MenuAction("Show Header Panel", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                if (menuItem.isSelected()) {
                    IGV.getInstance().getMainPanel().restoreHeader();
                } else {
                    IGV.getInstance().getMainPanel().removeHeader();
                }
                IGV.getInstance().repaint();
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction, true));

        menuItems.add(new JSeparator());
        menuAction =
                new MenuAction("Reorder Panels...", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        ReorderPanelsDialog dlg = new ReorderPanelsDialog(IGV.getMainFrame());
                        dlg.setVisible(true);
                    }
                };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction =
                new MenuAction("Add New Panel", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        String newPanelName = "Panel" + System.currentTimeMillis();
                        IGV.getInstance().addDataPanel(newPanelName);
                    }
                };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());
        menuItems.add(new HistoryMenu("Go to"));


        // Add to IGVPanel menu
        MenuAction dataMenuAction = new MenuAction("View", null, KeyEvent.VK_V);
        viewMenu = MenuAndToolbarUtils.createMenu(menuItems, dataMenuAction);
        return viewMenu;
    }

    private JMenu createRegionsMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;


        menuAction = new NavigateRegionsMenuAction("Region Navigator ...", IGV.getInstance());
        menuAction.setToolTipText(UIConstants.REGION_NAVIGATOR_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction =
                new MenuAction("Gene Lists...", null, KeyEvent.VK_S) {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        (GeneListManagerUI.getInstance(IGV.getMainFrame())).setVisible(true);
                    }
                };
        menuAction.setToolTipText("Open gene list manager");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        // Export Regions
        menuAction = new ExportRegionsMenuAction("Export Regions ...", KeyEvent.VK_E, IGV.getInstance());
        menuAction.setToolTipText(UIConstants.EXPORT_REGION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        // Import Regions
        menuAction = new ImportRegionsMenuAction("Import Regions ...", KeyEvent.VK_I, IGV.getInstance());
        menuAction.setToolTipText(IMPORT_REGION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Import Regions
//         menuAction = new ClearRegionsMenuAction("Clear Regions ...", IGV.getInstance());
//         menuAction.setToolTipText(IMPORT_REGION_TOOLTIP);
//         menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        MenuAction dataMenuAction = new MenuAction("Regions", null, KeyEvent.VK_V);
        viewMenu = MenuAndToolbarUtils.createMenu(menuItems, dataMenuAction);
        return viewMenu;
    }

    private JMenu createHelpMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();

        MenuAction menuAction = null;

        menuAction =
                new MenuAction("User Guide ... ") {

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


        if (Desktop.isDesktopSupported()) {
            final Desktop desktop = Desktop.getDesktop();
            if (desktop.isSupported(Desktop.Action.MAIL)) {

                menuAction =
                        new MenuAction("Help Forum...") {

                            @Override
                            public void actionPerformed(ActionEvent e) {
                                try {
                                    URI uri = new URI("http://groups.google.com/forum/#!forum/igv-help");
                                    Desktop.getDesktop().browse(uri);
                                } catch (Exception ex) {
                                    log.error("Error opening igv-help uri", ex);
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
                        (new AboutDialog(IGV.getMainFrame(), true)).setVisible(true);
                    }
                };
        menuAction.setToolTipText(ABOUT_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        MenuAction helpMenuAction = new MenuAction("Help");


        return MenuAndToolbarUtils.createMenu(menuItems, helpMenuAction);
    }

    private void checkVersion() {

        int readTimeout = Globals.READ_TIMEOUT;
        int connectTimeout = Globals.CONNECT_TIMEOUT;

        try {
            Main.Version thisVersion = Main.Version.getVersion(Globals.VERSION);

            if (thisVersion != null) {

                Globals.CONNECT_TIMEOUT = 5000;
                Globals.READ_TIMEOUT = 1000;
                final String serverVersionString = HttpUtils.getInstance().getContentsAsString(HttpUtils.createURL(Globals.getVersionURL())).trim();
                // See if user has specified to skip this update

                final String skipString = PreferencesManager.getPreferences().get(SKIP_VERSION);
                HashSet<String> skipVersion = new HashSet<String>(Arrays.asList(skipString.split(",")));
                if (skipVersion.contains(serverVersionString)) return;

                Main.Version serverVersion = Main.Version.getVersion(serverVersionString.trim());
                if (serverVersion == null) return;

                if (thisVersion.lessThan(serverVersion)) {

                    log.info("A later version of IGV is available (" + serverVersionString + ")");
                    final VersionUpdateDialog dlg = new VersionUpdateDialog(serverVersionString);

                    dlg.setVisible(true);
                    if (dlg.isSkipVersion()) {
                        String newSkipString = skipString + "," + serverVersionString;
                        PreferencesManager.getPreferences().put(SKIP_VERSION, newSkipString);
                    }

                } else {
                    MessageUtils.showMessage("IGV is up to date");
                }
            } else {
                if (Globals.VERSION.contains("3.0_beta") || Globals.VERSION.contains("snapshot")) {
                    HttpUtils.getInstance().getContentsAsString(HttpUtils.createURL(Globals.getVersionURL())).trim();
                }
            }

        } catch (Exception e) {
            log.error("Error checking version", e);
        } finally {
            Globals.CONNECT_TIMEOUT = connectTimeout;
            Globals.READ_TIMEOUT = readTimeout;
        }
    }

    private JMenu createExtrasMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();

        JMenuItem memTest = new JMenuItem("Memory test");
        memTest.addActionListener(e -> {
            CommandExecutor exe = new CommandExecutor(this.igv);
            int count = 1;
            int start = 0;
            exe.execute("snapshotDirectory /Users/jrobinso/Downloads/tmp");
            while (count++ < 10000) {
                exe.execute("goto chr1:" + start + "-" + (start + 1000));
                exe.execute("snapshot");
                start += 1000;
            }
        });
        menuItems.add(memTest);

        MenuAction menuAction = null;

        // Preferences reset
        menuAction = new ResetPreferencesAction("Reset Preferences", IGV.getInstance());
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
        menuItems.add(new JSeparator());


        // Set frame dimensions
        menuAction =
                new MenuAction("Set window dimensions", null, KeyEvent.VK_C) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        String value = JOptionPane.showInputDialog("Enter dimensions, e.g. 800x400");
                        if (value != null) {
                            String[] vals = value.split("x");
                            if (vals.length == 2) {
                                int w = Integer.parseInt(vals[0]);
                                int h = Integer.parseInt(vals[1]);
                                IGV.getMainFrame().setSize(w, h);
                            }
                        }
                    }
                };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Save entire window
        menuAction =
                new MenuAction("Save Screenshot ...", null, KeyEvent.VK_A) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        IGV.getInstance().saveImage(IGV.getInstance().getContentPane());

                    }
                };

        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menuAction = new ExportTrackNamesMenuAction("Export track names...", IGV.getInstance());
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menuAction = new MenuAction("Scatter Plot ...") {
            @Override
            public void actionPerformed(ActionEvent e) {
                final ReferenceFrame defaultFrame = FrameManager.getDefaultFrame();
                String chr = defaultFrame.getChrName();
                int start = (int) defaultFrame.getOrigin();
                int end = (int) defaultFrame.getEnd();
                int zoom = defaultFrame.getZoom();
                ScatterPlotUtils.openPlot(chr, start, end, zoom);
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        MenuAction extrasMenuAction = new MenuAction("Extras");
        JMenu menu = MenuAndToolbarUtils.createMenu(menuItems, extrasMenuAction);


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
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                            break;
                        }
                    }
                }
            });
            lfMenu.add(cb);
        }
        menu.add(lfMenu);

        menu.setVisible(false);


        return menu;
    }

    private JMenu createAWSMenu() throws IOException {

        JMenu menu = new JMenu("Amazon");

        final JMenuItem login = new JMenuItem("Login");
        login.addActionListener(e -> {
            try {
                OAuthProvider oauth = OAuthUtils.getInstance().getProvider("Amazon");
                oauth.openAuthorizationPage(); // should trigger and event and UI takes over
            } catch (Exception ex) {
                MessageUtils.showErrorMessage("Error fetching oAuth tokens.  See log for details", ex);
                log.error("Error fetching oAuth tokens", ex);
            }
        });
        login.setEnabled(false);
        menu.add(login);

        final JMenuItem logout = new JMenuItem("Logout");
        logout.addActionListener(e -> {
            OAuthProvider oauth = OAuthUtils.getInstance().getProvider("Amazon");
            oauth.logout();
        });
        logout.setEnabled(false);
        menu.add(logout);

        final JMenuItem loadS3 = new JMenuItem("Load from S3 bucket");
        loadS3.addActionListener(e -> {
            List<String> buckets = AmazonUtils.ListBucketsForUser();
            log.debug(buckets);

            UIUtilities.invokeOnEventThread(() -> {
                S3LoadDialog dlg = new S3LoadDialog(IGV.getMainFrame());
                dlg.setModal(true);
                dlg.setVisible(true);
                dlg.dispose();
            });
        });
        loadS3.setEnabled(false);
        menu.add(loadS3);

        menu.addMenuListener(new MenuListener() {
            @Override
            public void menuSelected(MenuEvent e) {
                Runnable runnable = () -> {
                    OAuthProvider oauth = OAuthUtils.getInstance().getProvider("Amazon");
                    boolean loggedIn = false;
                    loggedIn = oauth.isLoggedIn();
                    log.debug("MenuBar is user loggedIn?: " + loggedIn);

                    if (loggedIn) {
                        login.setText(oauth.getCurrentUserName());
                    } else {
                        login.setText("Login ...");
                    }
                    login.setEnabled(!loggedIn);
                    logout.setEnabled(loggedIn);
                    loadS3.setEnabled(loggedIn);
                };

                LongRunningTask.submit(runnable);
            }

            @Override
            public void menuDeselected(MenuEvent e) {

            }

            @Override
            public void menuCanceled(MenuEvent e) {

            }

        });

        return menu;
    }

    private JMenu createGoogleMenu() throws IOException {

        // Dynamically name menu - dwm08
        final OAuthProvider oauth = OAuthUtils.getInstance().getProvider();

        if (oauth != null) {  // TODO -- how do we know this is a google provider?
            oauth.setAuthProvider("Google");
            JMenu menu = new JMenu(oauth.getAuthProvider());

            final JMenuItem login = new JMenuItem("Login ... ");
            login.addActionListener(e -> {
                try {
                    oauth.openAuthorizationPage();
                } catch (Exception ex) {
                    MessageUtils.showErrorMessage("Error fetching oAuth tokens.  See log for details", ex);
                    log.error("Error fetching oAuth tokens", ex);
                }

            });
            //login.setEnabled(false);
            menu.add(login);


            final JMenuItem logout = new JMenuItem("Logout ");
            logout.addActionListener(e -> {
                oauth.logout();
                GoogleUtils.setProjectID(null);
            });
            logout.setEnabled(false);
            menu.add(logout);

            final JMenuItem projectID = new JMenuItem("Enter Project ID ...");
            projectID.addActionListener(e -> GoogleUtils.enterGoogleProjectID());
            menu.add(projectID);

            menu.addMenuListener(new MenuListener() {
                @Override
                public void menuSelected(MenuEvent e) {
                    Runnable runnable = () -> {
                        boolean loggedIn = OAuthUtils.getInstance().getProvider().isLoggedIn();

                        if (loggedIn) {
                            login.setText(oauth.getCurrentUserName());
                        } else {
                            login.setText("Login ...");
                        }
                        login.setEnabled(!loggedIn);
                        logout.setEnabled(loggedIn);
                    };

                    LongRunningTask.submit(runnable);
                }

                @Override
                public void menuDeselected(MenuEvent e) {

                }

                @Override
                public void menuCanceled(MenuEvent e) {

                }

            });
            return menu;
        } else {
            return null;
        }
    }

//    public void enableRemoveGenomes() {
//        if (removeImportedGenomeAction != null) {
//            removeImportedGenomeAction.setEnabled(true);
//        }
//    }

    public void resetSessionActions() {
        if (filterTracksAction != null) {
            filterTracksAction.resetTrackFilter();
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

    public JMenu getViewMenu() {
        return viewMenu;
    }

    final public void doExitApplication() {

        try {
            IGV.getInstance().saveStateForExit();
            Frame mainFrame = IGV.getMainFrame();
            PreferencesManager.getPreferences().setApplicationFrameBounds(mainFrame.getBounds());

            // Hide and close the application
            mainFrame.setVisible(false);
            mainFrame.dispose();

        } finally {
            System.exit(0);
        }
    }

    @ForTesting
    static void destroyInstance() {
        instance = null;
    }

    public void enableGoogleMenu(boolean aBoolean) {
        if (googleMenu != null) {
            googleMenu.setVisible(aBoolean);
        }
    }

    @Override
    public void receiveEvent(final Object event) {

        if (event instanceof GenomeChangeEvent) {
            UIUtilities.invokeOnEventThread(() -> encodeMenuItem.setVisible(EncodeFileBrowser.genomeSupported(((GenomeChangeEvent) event).genome.getId())));
        }
    }

    public void enableReloadSession() {
        this.reloadSessionItem.setEnabled(true);
    }

    public void disableReloadSession() {
        this.reloadSessionItem.setEnabled(false);
    }
}
