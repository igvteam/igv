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

import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.aws.S3LoadDialog;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.charts.ScatterPlotUtils;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.IGVEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.ChromSizesUtils;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.commandbar.HostedGenomeSelectionDialog;
import org.broad.igv.util.GoogleUtils;
import org.broad.igv.oauth.OAuthProvider;
import org.broad.igv.oauth.OAuthUtils;
import org.broad.igv.lists.GeneListManagerUI;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.tools.IgvToolsGui;
import org.broad.igv.tools.motiffinder.MotifFinderPlugin;
import org.broad.igv.track.CombinedDataSourceDialog;
import org.broad.igv.ui.action.*;
import org.broad.igv.ui.commandbar.RemoveGenomesDialog;
import org.broad.igv.ui.legend.LegendDialog;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.ReorderPanelsDialog;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.encode.EncodeTrackChooser;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import javax.swing.plaf.basic.BasicBorders;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.*;
import java.net.URI;
import java.util.*;
import java.util.List;
import java.util.function.Function;

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

    private static Logger log = LogManager.getLogger(IGVMenuBar.class);
    private static final String LOAD_GENOME_SERVER_TOOLTIP = "Select genomes available on the server to appear in menu.";
    private static final String CANNOT_LOAD_GENOME_SERVER_TOOLTIP = "Could not reach genome server";

    private static IGVMenuBar instance;

    private JMenu extrasMenu;
    private JMenu toolsMenu;
    private JMenu googleMenu;
    private JMenu AWSMenu;
    private AutosaveMenu autosaveMenu;
    private FilterTracksMenuAction filterTracksAction;
    private JMenu viewMenu;
    private IGV igv;
    /**
     * We store this as a field because we alter it if
     * we can't access genome server list
     */
    private JMenuItem loadGenomeFromServerMenuItem;
    private JMenuItem loadTracksFromServerMenuItem;
    private JMenuItem selectGenomeAnnotationsItem;

    private JMenuItem encodeUCSCMenuItem;
    private List<JComponent> encodeMenuItems = new ArrayList<>();

    private JMenuItem reloadSessionItem;
    private JMenuItem recentFilesMenu;


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


    public void notifyGenomeServerReachable(boolean reachable) {
        if (loadGenomeFromServerMenuItem != null) {
            UIUtilities.invokeOnEventThread(() -> {
                loadGenomeFromServerMenuItem.setEnabled(reachable);
                String tooltip = reachable ? LOAD_GENOME_SERVER_TOOLTIP : CANNOT_LOAD_GENOME_SERVER_TOOLTIP;
                loadGenomeFromServerMenuItem.setToolTipText(tooltip);
            });
        }
    }

    public void showAboutDialog() {
        (new AboutDialog(igv.getMainFrame(), true)).setVisible(true);
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

        // Create a placehold Google menu.  If not explicitly enabled it will remain invisible until triggered
        // by loading a protected Google resource
        try {
            googleMenu = createGoogleMenu();
            if (googleMenu != null) {
                boolean enabled = PreferencesManager.getPreferences().getAsBoolean(ENABLE_GOOGLE_MENU);
                enableGoogleMenu(enabled);
                menus.add(googleMenu);
            }
        } catch (IOException e) {
            log.error("Error creating google menu: " + e.getMessage());
        }


        AWSMenu = createAWSMenu();
        AWSMenu.setVisible(false);
        menus.add(AWSMenu);
        //detecting the provider is slow, do it in another thread
        LongRunningTask.submit(this::updateAWSMenu);


        menus.add(createHelpMenu());

        // Experimental -- remove for production release

        return menus;
    }

    public void updateAWSMenu() {
        SwingUtilities.invokeLater(() -> AWSMenu.setVisible(AmazonUtils.isAwsProviderPresent()));
    }

    /**
     * Generate the "tools" menu.
     * Legacy pattern -- at one times tools could be loaded dynamically as plug-ins
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
                IgvToolsGui.launch(false, GenomeManager.getInstance().getGenomeId());
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Motif finder
        menuItems.add(MotifFinderPlugin.getMenuItem());

        // BLAT
        menuItems.add(createBlatMenuItem());

        // Combine data tracks
        JMenuItem combineDataItem = new JMenuItem("Combine Data Tracks");
        combineDataItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                CombinedDataSourceDialog dialog = new CombinedDataSourceDialog(igv.getMainFrame());
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

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        String genomeId = genome == null ? null : genome.getId();

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

        menuAction = new LoadFromServerAction("Load From Server...", KeyEvent.VK_S, igv);
        menuAction.setToolTipText(UIConstants.LOAD_SERVER_DATA_TOOLTIP);
        loadTracksFromServerMenuItem = MenuAndToolbarUtils.createMenuItem(menuAction);
        menuItems.add(loadTracksFromServerMenuItem);

        recentFilesMenu = new RecentUrlsMenu();
        menuItems.add(recentFilesMenu);

        if (PreferencesManager.getPreferences().getAsBoolean(DB_ENABLED)) {
            menuAction = new LoadFromDatabaseAction("Load from Database...", 0, igv);
            menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
        }

        // ENCODE items.  These will be hidden / shown depending on genome chosen
        JSeparator separator = new JSeparator();
        menuItems.add(separator);

        JLabel encodeLabel = new JLabel("   ENCODE");
        encodeLabel.setFont(encodeLabel.getFont().deriveFont(Font.BOLD));
        menuItems.add(encodeLabel);

        // Post 2012 ENCODE menu
        JMenuItem chipItem = new JMenuItem();
        chipItem.setAction(new BrowseEncodeAction("Signals ChIP ...", 0, BrowseEncodeAction.Type.SIGNALS_CHIP, igv));
        encodeMenuItems.add(chipItem);

        JMenuItem otherSignalsItem = new JMenuItem();
        otherSignalsItem.setAction(new BrowseEncodeAction("Signals Other ...", 0, BrowseEncodeAction.Type.SIGNALS_OTHER, igv));
        encodeMenuItems.add(otherSignalsItem);

        JMenuItem otherItem = new JMenuItem();
        otherItem.setAction(new BrowseEncodeAction("Other ...", 0, BrowseEncodeAction.Type.OTHER, igv));
        encodeMenuItems.add(otherItem);

        for(JComponent item : encodeMenuItems) {
            menuItems.add(item);
            item.setVisible(EncodeTrackChooser.genomeSupported(genomeId));
        }

        // UCSC hosted ENCODE menu.
        encodeUCSCMenuItem = MenuAndToolbarUtils.createMenuItem(
                new BrowseEncodeAction("UCSC Repository (2012) ...", KeyEvent.VK_E, BrowseEncodeAction.Type.UCSC, igv));
        encodeUCSCMenuItem.setVisible(EncodeTrackChooser.genomeSupportedUCSC(genomeId));
        menuItems.add(encodeUCSCMenuItem);

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

        autosaveMenu = new AutosaveMenu();
        menuItems.add(autosaveMenu);

        menuItems.add(new JSeparator());

        // ***** Snapshots
        // Snapshot Application
        menuAction =
                new MenuAction("Save PNG Image ...", null, KeyEvent.VK_A) {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        igv.saveImage(igv.getMainPanel(), "png");

                    }
                };

        menuAction.setToolTipText(SAVE_PNG_IMAGE_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction =
                new MenuAction("Save SVG Image ...", null) {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        igv.saveImage(igv.getMainPanel(), "svg");

                    }
                };

        menuAction.setToolTipText(SAVE_SVG_IMAGE_TOOLTIP);
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
        JSeparator recentSessionsSep = new JSeparator();
        recentSessionsSep.setVisible(false);
        menuItems.add(recentSessionsSep);
        //menuItems.addAll(addRecentSessionMenuItems());
        menuItems.add(new JSeparator());;
        MenuAction fileMenuAction = new MenuAction("File", null, KeyEvent.VK_F);
        JMenu fileMenu = MenuAndToolbarUtils.createMenu(menuItems, fileMenuAction);


        //Add dynamic list of recent sessions
        fileMenu.addMenuListener(new DynamicItemsDisplay<>(fileMenu, recentSessionsSep, IGV.getInstance().getRecentSessionList(),
                session -> MenuAndToolbarUtils.createMenuItem(new OpenSessionMenuAction(session, IGV.getInstance()))));

        return fileMenu;
    }

    private void addEncodeItems(List<JComponent> menuItems, String genomeId) {

        JSeparator separator = new JSeparator();
        menuItems.add(separator);

        JLabel encodeLabel = new JLabel("   ENCODE");
        encodeLabel.setFont(encodeLabel.getFont().deriveFont(Font.BOLD));
        menuItems.add(encodeLabel);

        // Post 2012 ENCODE menu
        JMenuItem chipItem = new JMenuItem();
        chipItem.setAction(new BrowseEncodeAction("CHiP - Signals", 0, BrowseEncodeAction.Type.SIGNALS_CHIP, igv));
        encodeMenuItems.add(chipItem);

        JMenuItem otherSignalsItem = new JMenuItem();
        otherSignalsItem.setAction(new BrowseEncodeAction("Other - Signals", 0, BrowseEncodeAction.Type.SIGNALS_OTHER, igv));
        encodeMenuItems.add(otherSignalsItem);

        JMenuItem otherItem = new JMenuItem();
        otherItem.setAction(new BrowseEncodeAction("Other (peaks, calls, ...)", 0, BrowseEncodeAction.Type.OTHER, igv));
        encodeMenuItems.add(otherItem);

        for (JComponent item : encodeMenuItems) {
            menuItems.add(item);
            item.setVisible(EncodeTrackChooser.genomeSupported(genomeId));
        }

        // UCSC hosted ENCODE menu.
        encodeUCSCMenuItem = MenuAndToolbarUtils.createMenuItem(
                new BrowseEncodeAction("ENCODE (2012)...", KeyEvent.VK_E, BrowseEncodeAction.Type.UCSC, igv));
        encodeUCSCMenuItem.setVisible(EncodeTrackChooser.genomeSupportedUCSC(genomeId));
        menuItems.add(encodeUCSCMenuItem);
    }

    private JMenu createGenomesMenu() {

        JMenu menu = new JMenu("Genomes");

        loadGenomeFromServerMenuItem = new JMenuItem("Select Hosted Genome...");
        loadGenomeFromServerMenuItem.addActionListener(e -> HostedGenomeSelectionDialog.selectHostedGenome());
        loadGenomeFromServerMenuItem.setToolTipText(LOAD_GENOME_SERVER_TOOLTIP);
        menu.add(loadGenomeFromServerMenuItem);


        // Load genome json file
        JMenuItem fileItem = new JMenuItem("Load Genome from File...", KeyEvent.VK_I);
        fileItem.addActionListener(e -> {
            try {
                File importDirectory = PreferencesManager.getPreferences().getLastGenomeImportDirectory();
                if (importDirectory == null) {
                    PreferencesManager.getPreferences().setLastGenomeImportDirectory(DirectoryManager.getUserDirectory());
                }
                // Display the dialog
                File file = FileDialogUtils.chooseFile("Load Genome", importDirectory, FileDialog.LOAD);

                // If a file selection was made
                if (file != null) {
                    GenomeManager.getInstance().loadGenome(file.getAbsolutePath());
                }
            } catch (Exception ex) {
                MessageUtils.showErrorMessage(ex.getMessage(), ex);
            }
        });

        fileItem.setToolTipText("Load a FASTA, .json, or .genome file...");
        menu.add(fileItem);

        // Load genome from URL
        MenuAction urlMenuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_GENOME_FROM_URL, 0, igv);
        urlMenuAction.setToolTipText("Load a FASTA, .json, or .genome file...");
        menu.add(MenuAndToolbarUtils.createMenuItem(urlMenuAction));


        // Track hubs
        menu.add(new JSeparator());
        MenuAction genArkAction = new UCSCGenArkAction("Load Genome from UCSC GenArk...", 0, igv);
        menu.add(MenuAndToolbarUtils.createMenuItem(genArkAction));

        MenuAction menuAction = new SelectGenomeAnnotationTracksAction("Select GenArk Tracks...", igv);
        selectGenomeAnnotationsItem = MenuAndToolbarUtils.createMenuItem(menuAction);
        Genome newGenome = GenomeManager.getInstance().getCurrentGenome();
        selectGenomeAnnotationsItem.setEnabled(newGenome != null && newGenome.getHub() != null);
        menu.add(selectGenomeAnnotationsItem);
        menu.add(new JSeparator());

        // Add genome to combo box from server
        menuAction = new MenuAction("Remove Genomes...", null) {
            @Override
            public void actionPerformed(ActionEvent event) {
                RemoveGenomesDialog dialog2 = new RemoveGenomesDialog(igv.getMainFrame());
                dialog2.setVisible(true);
            }
        };
        menuAction.setToolTipText("Remove genomes which appear in the dropdown list");
        menu.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menu.addMenuListener(new MenuListener() {
            @Override
            public void menuSelected(MenuEvent e) {
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                selectGenomeAnnotationsItem.setEnabled(genome != null && genome.getHub() != null);
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

        // Export track names and attributes -- if > 1 i sselected export those, otherwise export all
        JMenuItem exportNames = new JMenuItem("Export Track Names and Attributes...");
        exportNames.addActionListener(e12 -> {
            Collection<Track> exportTracks = IGV.getInstance().getSelectedTracks();
            if (exportTracks.size() <= 1) {
                exportTracks = IGV.getInstance().getAllTracks();
            }
            exportTrackNames(exportTracks);
        });
        menuItems.add(exportNames);

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
                                igv.doViewPreferences();
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
                        (new LegendDialog(igv.getMainFrame())).setVisible(true);
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
                    igv.getMainPanel().expandNamePanel();
                } else {
                    igv.getMainPanel().collapseNamePanel();
                }
            }
        };
        boolean isShowing = igv.getMainPanel().isExpanded();
        JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
        menuItem.setSelected(isShowing);
        menuItem.setAction(menuAction);
        menuItems.add(menuItem);

        JMenuItem panelWidthmenuItem = new JMenuItem();
        menuAction = new MenuAction("Set Name Panel Width...", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {
                MainPanel mainPanel = igv.getMainPanel();
                String currentValue = String.valueOf(mainPanel.getNamePanelWidth());
                String newValue = MessageUtils.showInputDialog("Enter track name panel width: ", currentValue);
                if (newValue != null) {
                    try {
                        Integer w = Integer.parseInt(newValue);
                        if (w <= 0) throw new NumberFormatException();
                        PreferencesManager.getPreferences().put(NAME_PANEL_WIDTH, newValue);
                        mainPanel.setNamePanelWidth(w);
                    } catch (NumberFormatException ex) {
                        MessageUtils.showErrorMessage("Error: value must be a positive integer.", ex);
                    }
                }
            }
        };
        panelWidthmenuItem.setAction(menuAction);
        menuItems.add(panelWidthmenuItem);

        // Hide or Show the attribute panels
        //boolean isShow = PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY);
        //igv.doShowAttributeDisplay(isShow);

        menuAction = new MenuAction("Show Attribute Display", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                PreferencesManager.getPreferences().setShowAttributeView(menuItem.getState());
                igv.revalidateTrackPanels();
            }
        };
        boolean isShow = PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY);
        menuItem = MenuAndToolbarUtils.createMenuItem(menuAction, isShow);
        menuItems.add(menuItem);


        menuAction =
                new MenuAction("Select Attributes to Show...", null, KeyEvent.VK_S) {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        igv.doSelectDisplayableAttribute();
                    }
                };
        menuAction.setToolTipText(SELECT_DISPLAYABLE_ATTRIBUTES_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new MenuAction("Show Header Panel", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                if (menuItem.isSelected()) {
                    igv.getMainPanel().restoreHeader();
                } else {
                    igv.getMainPanel().removeHeader();
                }
                igv.getMainPanel().revalidate();
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction, true));

        menuItems.add(new JSeparator());
        menuAction =
                new MenuAction("Reorder Panels...", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        ReorderPanelsDialog dlg = new ReorderPanelsDialog(igv.getMainFrame());
                        dlg.setVisible(true);
                    }
                };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction =
                new MenuAction("Add New Panel", null, KeyEvent.VK_S) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        String newPanelName = "Panel" + System.currentTimeMillis();
                        igv.addDataPanel(newPanelName);
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
                        (GeneListManagerUI.getInstance(igv.getMainFrame())).setVisible(true);
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
                        (new AboutDialog(igv.getMainFrame(), true)).setVisible(true);
                    }
                };
        menuAction.setToolTipText(ABOUT_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        MenuAction helpMenuAction = new MenuAction("Help");


        return MenuAndToolbarUtils.createMenu(menuItems, helpMenuAction);
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
                                igv.getMainFrame().setSize(w, h);
                            }
                        }
                    }
                };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Save entire window
        menuAction =
                new MenuAction("Save PNG Screenshot ...", null, KeyEvent.VK_A) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        igv.saveImage(igv.getContentPane(), "png");

                    }
                };

        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction =
                new MenuAction("Save SVG Screenshot ...", null) {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        igv.saveImage(igv.getContentPane(), "svg");

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

        JMenuItem updateCS = new JMenuItem("Update chrom sizes");
        updateCS.addActionListener(e -> {
            try {
                ChromSizesUtils.main(new String[]{});
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
        menu.add(updateCS);

        menu.setVisible(false);
        return menu;
    }

    private JMenu createAWSMenu() {

        boolean usingCognito = AmazonUtils.GetCognitoConfig() != null;

        JMenu menu = new JMenu("Amazon");

        // Login
        final JMenuItem login = new JMenuItem("Login");
        login.addActionListener(e -> {
            try {
                OAuthProvider oauth = OAuthUtils.getInstance().getAWSProvider();
                oauth.openAuthorizationPage();
            } catch (Exception ex) {
                MessageUtils.showErrorMessage("Error fetching oAuth tokens.  See log for details", ex);
                log.error("Error fetching oAuth tokens", ex);
            }
        });
        login.setEnabled(usingCognito);
        login.setVisible(usingCognito);
        menu.add(login);

        // Logout
        final JMenuItem logout = new JMenuItem("Logout");
        logout.addActionListener(e -> {
            OAuthProvider oauth = OAuthUtils.getInstance().getAWSProvider();
            oauth.logout();
        });
        logout.setEnabled(false);
        logout.setVisible(usingCognito);
        menu.add(logout);

        // Load item, added to menu later
        final JMenuItem loadS3 = new JMenuItem("Load from S3 bucket");
        loadS3.addActionListener(e -> {
            List<String> buckets = AmazonUtils.ListBucketsForUser();
            log.debug(buckets);

            UIUtilities.invokeOnEventThread(() -> {
                S3LoadDialog dlg = new S3LoadDialog(igv.getMainFrame());
                dlg.setModal(true);
                dlg.setVisible(true);
                dlg.dispose();
            });
        });
        loadS3.setEnabled(!usingCognito);  // If using Cognito, disalbe initially
        menu.add(loadS3);

        menu.addMenuListener(new MenuListener() {
            @Override
            public void menuSelected(MenuEvent e) {
                if (AmazonUtils.GetCognitoConfig() != null) {
                    Runnable runnable = () -> {
                        OAuthProvider oauth = OAuthUtils.getInstance().getAWSProvider();
                        boolean loggedIn = oauth.isLoggedIn();
                        log.debug("MenuBar is user loggedIn?: " + loggedIn);

                        if (loggedIn) {
                            login.setText(oauth.getCurrentUserName());
                        } else {
                            login.setText("Login ...");
                        }
                        login.setVisible(true);
                        logout.setVisible(true);
                        login.setEnabled(!loggedIn);
                        logout.setEnabled(loggedIn);
                        loadS3.setEnabled(loggedIn);
                    };
                    LongRunningTask.submit(runnable);
                }
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

    private JMenu createGoogleMenu() {

        final OAuthProvider googleProvider = OAuthUtils.getInstance().getGoogleProvider();
        if (googleProvider == null) {
            log.error("Error creating google oauth provider");
            return null;
        }


        googleMenu = new JMenu("Google");

        final JMenuItem login = new JMenuItem("Login ... ");

        login.addActionListener(e -> {
            try {
                googleProvider.openAuthorizationPage();
            } catch (Exception ex) {
                MessageUtils.showErrorMessage("Error fetching oAuth tokens.  See log for details", ex);
                log.error("Error fetching oAuth tokens", ex);
            }

        });
        googleMenu.add(login);

        final JMenuItem logout = new JMenuItem("Logout ");
        logout.addActionListener(e -> {
            googleProvider.logout();
            GoogleUtils.setProjectID(null);
        });
        googleMenu.add(logout);

        final JMenuItem projectID = new JMenuItem("Enter Project ID ...");
        projectID.addActionListener(e -> GoogleUtils.enterGoogleProjectID());
        googleMenu.add(projectID);

        googleMenu.addMenuListener(new MenuListener() {
            @Override
            public void menuSelected(MenuEvent e) {
                boolean loggedIn = googleProvider.isLoggedIn();
                if (loggedIn && googleProvider.getCurrentUserName() != null) {
                    login.setText(googleProvider.getCurrentUserName());
                } else {
                    login.setText("Login ...");
                }
                login.setEnabled(!loggedIn);
                logout.setEnabled(loggedIn);
            }

            @Override
            public void menuDeselected(MenuEvent e) {

            }

            @Override
            public void menuCanceled(MenuEvent e) {

            }

        });

        return googleMenu;
    }


    /**
     * The Google menu is enabled dynamically to defer loading of oAuth properties until needed.
     * *
     *
     * @return
     * @throws IOException
     */
    public void enableGoogleMenu(boolean enable) throws IOException {
        if (googleMenu != null) {
            googleMenu.setVisible(enable);
        }
    }

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

    final public void doExitApplication() {

        try {
            igv.saveStateForExit();
            Frame mainFrame = igv.getMainFrame();
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


    @Override
    public void receiveEvent(final IGVEvent event) {

        if (event instanceof GenomeChangeEvent) {
            UIUtilities.invokeOnEventThread(() -> {
                final Genome genome = ((GenomeChangeEvent) event).genome();
                final String genomeId = genome.getId();
                encodeUCSCMenuItem.setVisible(EncodeTrackChooser.genomeSupportedUCSC(genomeId));
                for (JComponent item : encodeMenuItems) {
                    item.setVisible(EncodeTrackChooser.genomeSupported(genomeId));
                }

            });
        }
    }

    public void enableReloadSession() {
        this.reloadSessionItem.setEnabled(true);
    }

    public void showRecentFilesMenu(){
        this.recentFilesMenu.setVisible(true);
    }

    public void disableReloadSession() {
        this.reloadSessionItem.setEnabled(false);
    }

    public static JMenuItem createBlatMenuItem() {
        JMenuItem menuItem = new JMenuItem("BLAT ...");
        menuItem.addActionListener(e -> {

            String blatSequence = MessageUtils.showInputDialog("Enter sequence to blat:");
            if (blatSequence != null) {
                if (blatSequence.length() < 20 || blatSequence.length() > 8000) {
                    MessageUtils.showMessage("BLAT sequences must be between 20 and 8000 bases in length.");
                } else {
                    BlatClient.doBlatQuery(blatSequence, "BLAT");
                }
            }
        });

        return menuItem;
    }

    private void exportTrackNames(final Collection<Track> selectedTracks) {

        if (selectedTracks.isEmpty()) {
            return;
        }

        File file = FileDialogUtils.chooseFile("Export track names",
                PreferencesManager.getPreferences().getLastTrackDirectory(),
                new File("trackNames.txt"),
                FileDialogUtils.SAVE);

        if (file == null) {
            return;
        }

        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)))) {

            List<String> attributes = AttributeManager.getInstance().getVisibleAttributes();

            pw.print("Name");
            for (String att : attributes) {
                pw.print("\t" + att);
            }
            pw.println();

            for (Track track : selectedTracks) {
                //We preserve the alpha value. This is motivated by MergedTracks
                pw.print(track.getName());

                for (String att : attributes) {
                    String val = track.getAttributeValue(att);
                    pw.print("\t" + (val == null ? "" : val));
                }
                pw.println();
            }
        } catch (IOException e) {
            MessageUtils.showErrorMessage("Error writing to file", e);
            log.error(e);
        }
    }

    private static class DynamicItemsDisplay<T> implements MenuListener {
        private final Collection<T> values;
        private final Function<T, JMenuItem> itemConstructor;
        private final List<JComponent> activeComponents;
        private final JMenu menu;
        private final JSeparator recentSessionsSep;

        public DynamicItemsDisplay(JMenu menu, JSeparator insertionPoint, Collection<T> values, Function<T, JMenuItem> itemConstructor) {
            this.menu = menu;
            this.recentSessionsSep = insertionPoint;
            this.values = values;
            this.itemConstructor = itemConstructor;
            this.activeComponents = new ArrayList<>();
        }

        private List<JMenuItem> getCurrentItems() {
            return values.stream().map(itemConstructor).toList();
        }

        @Override
        public void menuSelected(MenuEvent e) {
            List<JMenuItem> newComponents = getCurrentItems();
            activeComponents.forEach(menu::remove);
            if (newComponents.isEmpty()) {
                recentSessionsSep.setVisible(false);
            } else {
                recentSessionsSep.setVisible(true);
                final int componentIndex = Arrays.asList(menu.getMenuComponents()).indexOf(recentSessionsSep);
                for(int i = 0; i < newComponents.size(); i++ ){
                    menu.insert(newComponents.get(i), componentIndex+i+1);
                }
                activeComponents.addAll(newComponents);
            }
            menu.revalidate();
            menu.repaint();
        }

        @Override
        public void menuDeselected(MenuEvent e) {}

        @Override
        public void menuCanceled(MenuEvent e) {}
    }
}

