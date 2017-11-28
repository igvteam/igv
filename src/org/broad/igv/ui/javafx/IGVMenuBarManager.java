/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2017 Broad Institute
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
package org.broad.igv.ui.javafx;

import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.control.*;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.apache.log4j.Logger;
import org.broad.igv.ui.javafx.panel.MainContentPane;
import org.broad.igv.util.FileUtils;

import java.io.File;

// Intended as the rough equivalent of the IGVMenuBar class of the Swing UI.  Work in progress.
// Will add event handlers (or at least stubs) for all of the included controls.
public class IGVMenuBarManager {
    private static Logger log = Logger.getLogger(IGVMenuBarManager.class);
    
    private MenuBar menuBar;

    // Keep as instance var for later break-out of actions, etc from constructor.
    private MainContentPane mainContentPane;

    public IGVMenuBarManager(Stage stage, MainContentPane mainContentPane) {
        this.mainContentPane = mainContentPane;

        // I'm leaving the creation of all of these inline for now. Need to break them
        // out for structural purposes and
        // to hold them as instance vars in order to manage enable/disable, handle
        // events, etc. We're not there yet so
        // it's too early to tell what is the best structure.

        // TODO: add actions to all of these MenuItems
        MenuItem loadFromFile = new MenuItem("Load from File ...");
        loadFromFile.setOnAction(new EventHandler<ActionEvent>() {

            @Override
            public void handle(ActionEvent actionEvent) {
                // Testing a FileChooser but with no real action
                FileChooser fileChooser = new FileChooser();
                fileChooser.setTitle("Choose a file");
                fileChooser.showOpenDialog(stage);
            }
        });
        MenuItem loadFromURL = new MenuItem("Load from URL ...");
        MenuItem loadFromServer = new MenuItem("Load from Server ...");
        MenuItem loadFromGa4gh = new MenuItem("Load from Ga4gh ...");
        MenuItem newSession = new MenuItem("New Session ...");
        MenuItem openSession = new MenuItem("Open Session ...");
        openSession.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent actionEvent) {
                // TODO: file filtering?
                FileChooser fileChooser = new FileChooser();
                fileChooser.setTitle("Choose a session file");
                File selected = fileChooser.showOpenDialog(stage);
                if (selected != null) {
                    String sessionFile = selected.getAbsolutePath();
                    log.info("About to load session");
                    if (sessionFile != null) {
                        if (FileUtils.isRemote(sessionFile)) {
                            // TODO: we are not really dealing with the Session yet, so this is a placeholder for now.
                            //boolean merge = false;
                            //IGV.getInstance().doRestoreSession(sessionFile, null, merge);
                        } else {
                            File f = new File(sessionFile);
                            //IGV.getInstance().doRestoreSession(f, null);
                        }
                    }
                    log.info("Session loading underway");
                }
            }
        });
        MenuItem saveSession = new MenuItem("Save Session ...");
        MenuItem saveImage = new MenuItem("Save Image ...");
        MenuItem exit = new MenuItem("Exit");
        exit.setOnAction(e -> Platform.exit());
        
        Menu fileMenu = new Menu("File", null, loadFromFile, loadFromURL, loadFromServer, loadFromGa4gh,
                new SeparatorMenuItem(), newSession, openSession, saveSession, new SeparatorMenuItem(), saveImage,
                new SeparatorMenuItem(), exit);

        MenuItem loadGenomeFromFile = new MenuItem("Load Genome From File ...");
        MenuItem loadGenomeFromURL = new MenuItem("Load Genome From URL ...");
        MenuItem loadGenomeFromServer = new MenuItem("Load Genome From Server ...");
        MenuItem createDotGenomeFile = new MenuItem("Create .genome File ...");
        MenuItem manageGenomeList = new MenuItem("Manage Genome List ...");
        Menu genomesMenu = new Menu("Genomes", null, loadGenomeFromFile, loadGenomeFromURL, loadGenomeFromServer,
                new SeparatorMenuItem(), createDotGenomeFile, new SeparatorMenuItem(), manageGenomeList);

        MenuItem preferences = new MenuItem("Preferences");
        MenuItem colorLegends = new MenuItem("Color Legends ...");
        CheckMenuItem showNamePanel = new CheckMenuItem("Show Name Panel");
        showNamePanel.setSelected(true);
        showNamePanel.setOnAction(new EventHandler<ActionEvent>() {

            @Override
            public void handle(ActionEvent event) {
                if (mainContentPane != null) {
                    if (mainContentPane.isNamePanelHidden()) {
                        mainContentPane.showNamePanel();
                    } else {
                        mainContentPane.hideNamePanel();
                    }
                }
            }
        });
        MenuItem setNamePanelWidth = new MenuItem("Set Name Panel Width ...");
        CheckMenuItem showAttribsDisplay = new CheckMenuItem("Show Attributes Display");
        MenuItem selectAttribsToShow = new MenuItem("Select Attributes to Show ...");
        CheckMenuItem showHeaderPanel = new CheckMenuItem("Show Header Panel");
        MenuItem reorderPanels = new MenuItem("Reorder Panels ...");

        MenuItem gotoBack = new MenuItem("Back");
        MenuItem gotoFwd = new MenuItem("Forward");
        MenuItem clearAll = new MenuItem("Clear All");
        Menu gotoSubMenu = new Menu("Go to", null, gotoBack, gotoFwd, new SeparatorMenuItem(), clearAll);

        Menu viewMenu = new Menu("View", null, preferences, colorLegends, new SeparatorMenuItem(), showNamePanel,
                setNamePanelWidth, showAttribsDisplay, selectAttribsToShow, showHeaderPanel, new SeparatorMenuItem(),
                reorderPanels, new SeparatorMenuItem(), gotoSubMenu, new SeparatorMenuItem());

        MenuItem sortTracks = new MenuItem("Sort Tracks ...");
        MenuItem groupTracks = new MenuItem("Group Tracks ...");
        MenuItem filterTracks = new MenuItem("Filter Tracks ...");
        MenuItem fitDataToWindow = new MenuItem("Fit Data to Window");
        MenuItem setTrackHeight = new MenuItem("Set Track Height ...");

        Menu tracksMenu = new Menu("Tracks", null, sortTracks, groupTracks, filterTracks, new SeparatorMenuItem(),
                fitDataToWindow, setTrackHeight);

        MenuItem navigateRegions = new MenuItem("Navigate Regions ...");
        MenuItem geneLists = new MenuItem("Gene Lists ...");
        MenuItem exportRegions = new MenuItem("Export Regions ...");
        MenuItem importRegions = new MenuItem("Import Regions ...");
        Menu regionsMenu = new Menu("Regions", null, navigateRegions, geneLists, new SeparatorMenuItem(), exportRegions,
                importRegions);

        MenuItem runBatchScript = new MenuItem("Run Batch Script ...");
        MenuItem runIgvtools = new MenuItem("Run igvtools ...");
        MenuItem findMotif = new MenuItem("Find Motif ...");
        MenuItem blat = new MenuItem("BLAT ...");
        MenuItem combineDataTracks = new MenuItem("Combine Data Tracks");

        MenuItem loadGeneMatrixInGitools = new MenuItem("Load GeneMatrix in Gitools");
        MenuItem exportGeneMatrix = new MenuItem("Export GeneMatrix (TDM) ...");
        Menu gitoolsHeatmaps = new Menu("Gitools Heatmaps", null, loadGeneMatrixInGitools, exportGeneMatrix);

        MenuItem intersect = new MenuItem("Intersect");
        MenuItem removeSubtract = new MenuItem("Remove/Subtract");
        MenuItem closest = new MenuItem("Closest");
        MenuItem window = new MenuItem("Windows");
        MenuItem coverage = new MenuItem("Coverage");
        MenuItem multiIntersect = new MenuItem("Multi-intersect");
        MenuItem setPathToBEDTools = new MenuItem("Set Path to BEDTools ...");
        Menu bedTools = new Menu("BEDTools", null, intersect, removeSubtract, closest, window, coverage, multiIntersect,
                setPathToBEDTools);
        Menu toolsMenu = new Menu("Tools", null, runBatchScript, runIgvtools, findMotif, blat, combineDataTracks,
                new SeparatorMenuItem(), gitoolsHeatmaps, bedTools);

        MenuItem loadFileFromGS = new MenuItem("Load File From GenomeSpace ...");
        MenuItem loadGenomeFromGS = new MenuItem("Load Genome From GenomeSpace ...");
        MenuItem saveSessionToGS = new MenuItem("Save Session To GenomeSpace ...");
        MenuItem loadSessionFromGS = new MenuItem("Load Session From GenomeSpace ...");
        MenuItem logout = new MenuItem("Logout");
        MenuItem register = new MenuItem("Register ...");

        Menu genomeSpaceMenu = new Menu("GenomeSpace", null, loadFileFromGS, new SeparatorMenuItem(), loadGenomeFromGS,
                new SeparatorMenuItem(), saveSessionToGS, loadSessionFromGS, new SeparatorMenuItem(), logout,
                new SeparatorMenuItem(), register);

        MenuItem userGuide = new MenuItem("User Guide ...");
        MenuItem helpForum = new MenuItem("Help Forum ...");
        MenuItem checkForUpdates = new MenuItem("Check For Updates ...");
        MenuItem aboutIGV = new MenuItem("About IGV");
        Menu helpMenu = new Menu("Help", null, userGuide, helpForum, checkForUpdates, aboutIGV);

        // TODO: need to add Google Menu, DB menu (?), Extras Menu, L&F (Skins) menu.
        // Maybe.

        menuBar = new MenuBar(fileMenu, genomesMenu, viewMenu, tracksMenu, regionsMenu, toolsMenu, genomeSpaceMenu,
                helpMenu);
        final String os = System.getProperty("os.name");
        if (os != null && os.startsWith("Mac"))
            menuBar.useSystemMenuBarProperty().set(true);
    }

    public MenuBar getMenuBar() {
        return menuBar;
    }
}
