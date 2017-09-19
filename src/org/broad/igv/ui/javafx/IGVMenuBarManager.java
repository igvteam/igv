package org.broad.igv.ui.javafx;

import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import javafx.stage.FileChooser;

import org.broad.igv.ui.util.UIUtilities;

// Intended as the rough equivalent of the IGVMenuBar class of the Swing UI.  Work in progress.
public class IGVMenuBarManager {
  private static IGVMenuBarManager instance;

  // TODO: decide if we need to capture this.  It's a singleton, so it's always available via its getInstance() method.
  IGVStageManager igv;

  private MenuBar menuBar;
  
  static IGVMenuBarManager createInstance(IGVStageManager igv) {
      if (instance != null) {
          if (igv == instance.igv) {
              return instance;
          }
          throw new IllegalStateException("Cannot create another IGVMenuBarManager, use getInstance");
      }
      UIUtilities.invokeAndWaitOnEventThread(() ->instance = new IGVMenuBarManager(igv));
      return instance;
  }

  public static IGVMenuBarManager getInstance() {
      return instance;
  }

  private IGVMenuBarManager(IGVStageManager igv) {
    // I'm leaving the creation of all of these inline for now.  Need to break them out for structural purposes and 
    // to hold them as instance vars in order to manage enable/disable, handle events, etc.  We're not there yet so
    // it's too early to tell what is the best structure.
    
    // TODO: add actions to all of these MenuItems
    MenuItem loadFromFile = new MenuItem("Load from File ...");
    loadFromFile.setOnAction(new EventHandler<ActionEvent>() {
      
      @Override
      public void handle(ActionEvent actionEvent) {
        // Testing a FileChooser but with no real action
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Choose a file");
        fileChooser.showOpenDialog(igv.getMainStage());
      }
    });
    MenuItem loadFromURL = new MenuItem("Load from URL ...");
    MenuItem loadFromServer = new MenuItem("Load from Server ...");
    MenuItem loadFromGa4gh = new MenuItem("Load from Ga4gh ...");
    MenuItem newSession = new MenuItem("New Session ...");
    MenuItem openSession = new MenuItem("Open Session ...");
    MenuItem saveSession = new MenuItem("Save Session ...");
    MenuItem saveImage = new MenuItem("Save Image ...");
    MenuItem exit = new MenuItem("Exit");
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
    MenuItem setNamePanelWidth = new MenuItem("Set Name Panel Width ...");
    CheckMenuItem showAttribsDisplay = new CheckMenuItem("Show Attributes Display");
    MenuItem selectAttribsToShow = new MenuItem("Select Attributes to Show ...");
    CheckMenuItem showHeaderPanel = new CheckMenuItem("Show Header Panel");
    MenuItem reorderPanels = new MenuItem("Reorder Panels ...");

    MenuItem gotoBack = new MenuItem("Back");
    MenuItem gotoFwd = new MenuItem("Forward");
    MenuItem clearAll = new MenuItem("Clear All");
    Menu gotoSubMenu = new Menu("Go to", null, gotoBack, gotoFwd, new SeparatorMenuItem(), clearAll);
    
    Menu viewMenu = new Menu("View", null, preferences, colorLegends, new SeparatorMenuItem(),
        showNamePanel, setNamePanelWidth, showAttribsDisplay, selectAttribsToShow, showHeaderPanel,
        new SeparatorMenuItem(), reorderPanels, new SeparatorMenuItem(), gotoSubMenu);
    
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
    Menu regionsMenu = new Menu("Regions", null, navigateRegions, geneLists, new SeparatorMenuItem(), exportRegions, importRegions);
    
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
    MenuItem multiIntersect = new MenuItem("Multi-interesect");
    MenuItem setPathToBEDTools = new MenuItem("Set Path to BEDTools ...");
    Menu bedTools = new Menu("BEDTools", null, intersect, removeSubtract, closest, window, coverage, multiIntersect, setPathToBEDTools);
    Menu toolsMenu = new Menu("Tools", null, runBatchScript, runIgvtools, findMotif, blat, combineDataTracks,
        new SeparatorMenuItem(), gitoolsHeatmaps, bedTools);


    MenuItem loadFileFromGS = new MenuItem("Load File From GenomeSpace ...");
    MenuItem loadGenomeFromGS = new MenuItem("Load Genome From GenomeSpace ...");
    MenuItem saveSessionToGS = new MenuItem("Save Session To GenomeSpace ...");
    MenuItem loadSessionFromGS = new MenuItem("Load Session From GenomeSpace ...");
    MenuItem logout = new MenuItem("Logout");
    MenuItem register = new MenuItem("Register ...");
    
    Menu genomeSpaceMenu = new Menu("GenomeSpace", null, loadFileFromGS, new SeparatorMenuItem(), loadGenomeFromGS, 
        new SeparatorMenuItem(), saveSessionToGS, loadSessionFromGS, new SeparatorMenuItem(), logout, new SeparatorMenuItem(), register);

    MenuItem userGuide = new MenuItem("User Guide ...");
    MenuItem helpForum = new MenuItem("Help Forum ...");
    MenuItem checkForUpdates = new MenuItem("Check For Updates ...");
    MenuItem aboutIGV = new MenuItem("About IGV");
    Menu helpMenu = new Menu("Help", null, userGuide, helpForum, checkForUpdates, aboutIGV);
    
    // TODO: need to add Google Menu, DB menu (?), Extras Menu, L&F (Skins) menu.  Maybe.
    
    menuBar = new MenuBar(fileMenu, genomesMenu, viewMenu, tracksMenu, regionsMenu, toolsMenu, genomeSpaceMenu, helpMenu);
    final String os = System.getProperty ("os.name");
    if (os != null && os.startsWith ("Mac"))
      menuBar.useSystemMenuBarProperty ().set (true);
  }

  public MenuBar getMenuBar() {
    return menuBar;
  }
}
