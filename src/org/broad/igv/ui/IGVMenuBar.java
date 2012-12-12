/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui;

import apple.dts.samplecode.osxadapter.OSXAdapter;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.charts.ScatterPlotUtils;
import org.broad.igv.dev.plugin.PluginSpecReader;
import org.broad.igv.dev.plugin.ui.RunPlugin;
import org.broad.igv.dev.plugin.ui.SetPluginPathDialog;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.gs.GSOpenSessionMenuAction;
import org.broad.igv.gs.GSSaveSessionMenuAction;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.lists.GeneListManagerUI;
import org.broad.igv.lists.VariantListManager;
import org.broad.igv.tools.IgvToolsGui;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.action.*;
import org.broad.igv.ui.legend.LegendDialog;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.ReorderPanelsDialog;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.BrowserLauncher;
import org.broad.tribble.Feature;
import org.w3c.dom.Element;

import javax.swing.*;
import javax.swing.plaf.basic.BasicBorders;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import static org.broad.igv.ui.UIConstants.*;

/**
 * @author jrobinso
 * @date Apr 4, 2011
 */
public class IGVMenuBar extends JMenuBar {

    private static Logger log = Logger.getLogger(IGVMenuBar.class);
    public static final String GENOMESPACE_REG_TOOLTIP = "Register for GenomeSpace";
    public static final String GENOMESPACE_REG_PAGE = "http://www.genomespace.org/register";

    private JMenu extrasMenu;
    //private RemoveUserDefinedGenomeMenuAction removeImportedGenomeAction;
    private FilterTracksMenuAction filterTracksAction;
    private JMenu viewMenu;
    IGV igv;

    private JMenu toolsMenu;

    public void showAboutDialog() {
        (new AboutDialog(IGV.getMainFrame(), true)).setVisible(true);
    }

    public IGVMenuBar(IGV igv) {
        this.igv = igv;
        setBorder(new BasicBorders.MenuBarBorder(Color.GRAY, Color.GRAY));
        setBorderPainted(true);

        for (AbstractButton menu : createMenus()) {
            add(menu);
        }

        //This is for Macs, so showing the about dialog
        //from the command bar does what we want.
        if (Globals.IS_MAC) {
            try {
                OSXAdapter.setAboutHandler(this, getClass().getDeclaredMethod("showAboutDialog", (Class[]) null));
                OSXAdapter.setQuitHandler(ShutdownThread.class, ShutdownThread.class.getDeclaredMethod("runS", (Class[]) null));
            } catch (Exception e) {
                log.error("Error setting apple-specific about and quit handlers", e);
            }

        }
    }

    private List<AbstractButton> createMenus() {

        List<AbstractButton> menus = new ArrayList<AbstractButton>();
        menus.add(createFileMenu());

        boolean affectiveMode = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.AFFECTIVE_ENABLE);
        if (!affectiveMode) {
            menus.add(createGenomesMenu());
        }

        menus.add(createViewMenu());
        menus.add(createTracksMenu());
        menus.add(createRegionsMenu());

        if (true || Globals.toolsMenuEnabled) {
            refreshToolsMenu();
            menus.add(toolsMenu);
        }

        menus.add(createGenomeSpaceMenu());
        extrasMenu = createExtrasMenu();
        //extrasMenu.setVisible(false);
        menus.add(extrasMenu);

        menus.add(createHelpMenu());

        // Experimental -- remove for production release

        return menus;
    }

    /**
     * Generate the "tools" menu.
     * This is imperative, it is written to field {@code toolsMenu}.
     * Reason being, when we add (TODO remove)
     * a new tool, we need to refresh just this menu
     */
    private void refreshToolsMenu() {
        List<JComponent> menuItems = new ArrayList<JComponent>(10);

        // batch script
        MenuAction menuAction = new RunScriptMenuAction("Run Batch Script...", KeyEvent.VK_X, igv);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // igvtools
        menuItems.add(new JSeparator());
        menuAction = new SortTracksMenuAction("Run igvtools...", KeyEvent.VK_T, igv) {
            @Override
            public void actionPerformed(ActionEvent e) {
                IgvToolsGui.launch(false, igv.getGenomeManager().getGenomeId());
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
        menuItems.add(new JSeparator());

        //-------------------------------------//
        //"Add tool" option, for loading plugin from someplace else
        JMenuItem addTool = new JMenuItem("Add tool");
        addTool.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                File pluginFi = FileDialogUtils.chooseFile("Select plugin .xml spec");
                if (pluginFi == null) return;

                try {
                    PluginSpecReader.addCustomPlugin(pluginFi.getAbsolutePath());
                    refreshToolsMenu();
                } catch (IOException e1) {
                    MessageUtils.showErrorMessage("Error loading custom plugin", e1);
                }
            }
        });
        //menuItems.add(addTool);
        //menuItems.add(new JSeparator());

        //-------------------------------------//

        JMenuItem exportData = new JMenuItem("Export Features");
        exportData.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                File outFile = FileDialogUtils.chooseFile("Save Visible Data",
                        PreferenceManager.getInstance().getLastTrackDirectory(),
                        new File("visibleData.bed"),
                        FileDialogUtils.SAVE);
                IGVMenuBar.exportVisibleData(outFile.getAbsolutePath(), IGV.getInstance().getAllTracks());
            }
        });
        //menuItems.add(exportData);

        //-------------------------------------//

        for (final PluginSpecReader pluginSpecReader : PluginSpecReader.getPlugins()) {
            for (final Element tool : pluginSpecReader.getTools()) {
                final String toolName = tool.getAttributes().getNamedItem("name").getTextContent();
                boolean toolVisible = Boolean.parseBoolean(tool.getAttribute("visible"));
                JMenuItem toolMenu;

                if (toolVisible) {

                    final String toolPath = pluginSpecReader.getToolPath(tool);
                    final String tool_url = tool.getAttribute("tool_url");
                    boolean isValid = PluginSpecReader.isToolPathValid(toolPath);

                    ActionListener invalidActionListener = new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            String msg = String.format("%s executable not found at %s", toolName, toolPath);
                            if (tool_url != null) {
                                msg += "<br/>See " + tool_url + " to install";
                            }
                            MessageUtils.showMessage(msg);
                        }
                    };

                    toolMenu = new JMenu(toolName);
                    //Kind of overlaps with the side-pull menu, doesn't look great
                    //toolMenu.setToolTipText(tool.getAttribute("description"));
                    for (final Element command : pluginSpecReader.getCommands(tool)) {
                        final String cmdName = command.getAttribute("name");
                        JMenuItem cmdItem = new JMenuItem(cmdName);
                        toolMenu.add(cmdItem);
                        if (isValid || toolPath.length() == 0) {
                            cmdItem.addActionListener(new ActionListener() {
                                @Override
                                public void actionPerformed(ActionEvent e) {
                                    (new RunPlugin(IGV.getMainFrame(), pluginSpecReader, tool, command)).setVisible(true);
                                }
                            });
                            cmdItem.setEnabled(true);
                        } else {
                            cmdItem.setEnabled(false);
                        }
                    }
                    if (toolPath.length() > 0) {
                        JMenuItem setPathItem = new JMenuItem(String.format("Set path to %s...", toolName));
                        setPathItem.addActionListener(new ActionListener() {
                            @Override
                            public void actionPerformed(ActionEvent e) {
                                (new SetPluginPathDialog(IGV.getMainFrame(), pluginSpecReader, tool)).setVisible(true);
                                refreshToolsMenu();
                            }
                        });
                        toolMenu.add(setPathItem);
                    }
                    menuItems.add(toolMenu);
                }
            }
        }

        //-------------------------------------//


        MenuAction toolsMenuAction = new MenuAction("Tools", null);
        if (toolsMenu == null) {
            toolsMenu = MenuAndToolbarUtils.createMenu(menuItems, toolsMenuAction);
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


    private JMenu createFileMenu() {

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

        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_FROM_DAS, KeyEvent.VK_D, igv);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.DB_ENABLED)) {
            menuAction = new LoadFromDatabaseAction("Load from Database...", 0, igv);
            menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));
        }

        menuItems.add(new JSeparator());

        // Session menu items
        menuAction = new NewSessionMenuAction("New Session...", KeyEvent.VK_N, igv);
        menuAction.setToolTipText(UIConstants.NEW_SESSION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new OpenSessionMenuAction("Open Session...", KeyEvent.VK_O, igv);
        menuAction.setToolTipText(UIConstants.RESTORE_SESSION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new SaveSessionMenuAction("Save Session...", KeyEvent.VK_V, igv);
        menuAction.setToolTipText(UIConstants.SAVE_SESSION_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

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
        String recentSessions = PreferenceManager.getInstance().getRecentSessions();
        if (recentSessions != null) {
            String[] sessions = recentSessions.split(";");
            for (String sessionPath : sessions) {
                if (!igv.getRecentSessionList().contains(sessionPath)) {
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
        return MenuAndToolbarUtils.createMenu(menuItems, fileMenuAction);
    }

    private void notifyGenomesAddedRemoved(List<GenomeListItem> selectedValues, boolean added) {
        if (selectedValues == null || selectedValues.size() == 0) return;
        int size = selectedValues.size();
        String msg = "";
        if (size == 1) {
            msg += selectedValues.get(0) + " genome";
        } else {
            msg += size + " genomes";
        }
        if (added) {
            msg += " added to";
        } else {
            msg += " removed from";
        }
        msg += " list";

        MessageUtils.setStatusBarMessage(msg);
    }


    private JMenu createGenomesMenu() {
        List<JComponent> menuItems = new ArrayList<JComponent>();
        MenuAction menuAction = null;

        // Load genome
        menuAction =
                new MenuAction("Load Genome from File...", null, KeyEvent.VK_I) {
                    @Override
                    public void actionPerformed(ActionEvent event) {
                        org.broad.igv.ui.util.ProgressMonitor monitor = new org.broad.igv.ui.util.ProgressMonitor();
                        igv.doLoadGenome(monitor);

                    }
                };

        menuAction.setToolTipText("Load a FASTA or .genome file...");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Load genome from URL
        menuAction = new LoadFromURLMenuAction(LoadFromURLMenuAction.LOAD_GENOME_FROM_URL, 0, igv);
        menuAction.setToolTipText("Load a FASTA or .genome file...");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        // Add genome to combo box from server
        menuAction = new MenuAction("Load Genome From Server...", null) {
            @Override
            public void actionPerformed(ActionEvent event) {
                GenomeSelectionDialog dialog = new GenomeSelectionDialog(IGV.getMainFrame(), ListSelectionModel.SINGLE_SELECTION);
                dialog.setVisible(true);
                List<GenomeListItem> selectedValues = dialog.getSelectedValuesList();
                if (selectedValues != null && selectedValues.size() >= 1) {
                    GenomeManager.getInstance().addGenomeItems(selectedValues);
                    igv.getContentPane().getCommandBar().refreshGenomeListComboBox();
                    //notifyGenomesAddedRemoved(selectedValues, true);
                    igv.selectGenomeFromList(selectedValues.get(0).getId());
                }
            }
        };
        menuAction.setToolTipText("Select genomes available on the server to appear in menu");
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        menuAction =
                new MenuAction("Create .genome File...", null, KeyEvent.VK_D) {
                    @Override
                    public void actionPerformed(ActionEvent event) {
                        org.broad.igv.ui.util.ProgressMonitor monitor = new org.broad.igv.ui.util.ProgressMonitor();
                        igv.doDefineGenome(monitor);
                    }
                };

        menuAction.setToolTipText(UIConstants.IMPORT_GENOME_TOOLTIP);
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        // Add genome to combo box from server
        menuAction = new MenuAction("Manage Genome List...", null) {
            @Override
            public void actionPerformed(ActionEvent event) {
                ManageGenomesDialog dialog2 = new ManageGenomesDialog(IGV.getMainFrame());
                dialog2.setVisible(true);
                boolean cancelled = dialog2.isCancelled();
                List<GenomeListItem> removedValuesList = dialog2.getRemovedValuesList();
                if (!cancelled) {
                    GenomeManager.getInstance().buildGenomeItemList();
                    igv.getContentPane().getCommandBar().refreshGenomeListComboBox();
                    if (removedValuesList != null && !removedValuesList.isEmpty()) {
                        GenomeManager.getInstance().updateImportedGenomePropertyFile();
                        notifyGenomesAddedRemoved(removedValuesList, false);
                    }
                }
            }
        };
        menuAction.setToolTipText("Add, remove, or reorder genomes which appear in the dropdown list");
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
                IGV.getInstance().doRefresh();
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
                        PreferenceManager.getInstance().put(PreferenceManager.NAME_PANEL_WIDTH, newValue);
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
        boolean isShow = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY);
        IGV.getInstance().doShowAttributeDisplay(isShow);  // <= WEIRD doing IGV.getInstance() here!

        menuAction = new MenuAction("Show Attribute Display", null, KeyEvent.VK_A) {
            @Override
            public void actionPerformed(ActionEvent e) {

                JCheckBoxMenuItem menuItem = (JCheckBoxMenuItem) e.getSource();
                PreferenceManager.getInstance().setShowAttributeView(menuItem.getState());
                IGV.getInstance().getMainPanel().invalidate();
                IGV.getInstance().doRefresh();
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
                IGV.getInstance().doRefresh();
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

        MenuAction toolMenuAction = new MenuAction("Help");
        return MenuAndToolbarUtils.createMenu(menuItems, toolMenuAction);
    }

    private JMenu createGenomeSpaceMenu() {

        JMenu menu = new JMenu("GenomeSpace");

        MenuAction menuAction = null;
        menuAction = new LoadFromGSMenuAction("Load from GenomeSpace...", KeyEvent.VK_U, IGV.getInstance());
        menu.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menu.addSeparator();

        menuAction = new GSSaveSessionMenuAction("Save session to GenomeSpace...", IGV.getInstance());
        menu.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuAction = new GSOpenSessionMenuAction("Load session from GenomeSpace...", IGV.getInstance());
        menu.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menu.add(new JSeparator());
        menuAction = new MenuAction("Logout") {
            @Override
            public void actionPerformed(ActionEvent e) {
                GSUtils.logout();
                if (MessageUtils.confirm("You must shutdown IGV to complete the GenomeSpace logout. Shutdown now?")) {
                    doExitApplication();
                }
            }
        };
        menu.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menu.add(new JSeparator());
        menuAction =
                new MenuAction("Register... ") {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        try {
                            BrowserLauncher.openURL(GENOMESPACE_REG_PAGE);
                        } catch (IOException ex) {
                            log.error("Error opening browser", ex);
                        }

                    }
                };
        menuAction.setToolTipText(GENOMESPACE_REG_TOOLTIP);
        menu.add(MenuAndToolbarUtils.createMenuItem(menuAction));


        menu.setVisible(PreferenceManager.getInstance().getAsBoolean(PreferenceManager.GENOME_SPACE_ENABLE));


        return menu;
    }

    private JMenu createExtrasMenu() {

        List<JComponent> menuItems = new ArrayList<JComponent>();

        MenuAction menuAction = null;

        final JCheckBoxMenuItem exomeModeItem = new JCheckBoxMenuItem("Exome mode");
        exomeModeItem.setSelected(FrameManager.isExomeMode());
        exomeModeItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                FrameManager.setExomeMode(exomeModeItem.isSelected(), true);
                igv.resetFrames();
            }
        });
        menuItems.add(exomeModeItem);
        menuItems.add(new JSeparator());


        // Preferences reset
        menuAction = new ResetPreferencesAction("Reset Preferences", IGV.getInstance());
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());

        menuAction = new MenuAction("Variant list ...  *EXPERIMENTAL*") {
            @Override
            public void actionPerformed(ActionEvent e) {
                VariantListManager.openNavigator(IGV.getMainFrame());
            }
        };
        menuItems.add(MenuAndToolbarUtils.createMenuItem(menuAction));

        menuItems.add(new JSeparator());


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
                            IGV.getMainFrame().setSize(w, h);
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
            // Hide and close the application
            mainFrame.setVisible(false);
            mainFrame.dispose();

        } finally {
            System.exit(0);
        }

    }

    /**
     * Write visible data to a file
     * TODO Move to own action class, thread
     */
    private static final void exportVisibleData(String outPath, Collection<Track> tracks) {
        PrintWriter writer;
        try {
            writer = new PrintWriter(outPath);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
        ReferenceFrame.Range range = FrameManager.getDefaultFrame().getCurrentRange();
        for (Track track : tracks) {
            if (track instanceof FeatureTrack) {
                FeatureTrack fTrack = (FeatureTrack) track;
                List<Feature> features = fTrack.getFeatures(range.getChr(), range.getStart(), range.getEnd());
                IGVBEDCodec codec = new IGVBEDCodec();
                for (Feature feat : features) {
                    String featString = codec.encode(feat);
                    writer.println(featString);
                }
            }
        }
        writer.flush();
        writer.close();
    }
}
