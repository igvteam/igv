package org.broad.igv.ui.javafx;

import javafx.application.Application;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Orientation;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Button;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.Slider;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TextField;
import javafx.scene.control.ToolBar;
import javafx.scene.image.ImageView;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.StrokeLineCap;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

public class PortingExperiment extends Application {

  public PortingExperiment() {
  }

  @Override
  public void start(Stage primaryStage) throws Exception {
    primaryStage.setTitle("JavaFX Porting Experiment");

    MenuItem loadFromFile = new MenuItem("Load from File ...");
    loadFromFile.setOnAction(new EventHandler<ActionEvent>() {
      
      @Override
      public void handle(ActionEvent arg0) {
        // Testing a FileChooser but with no real action
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Choose a file");
        fileChooser.showOpenDialog(primaryStage);
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
    
    MenuBar menuBar = new MenuBar(fileMenu, genomesMenu, viewMenu, tracksMenu, regionsMenu, toolsMenu, genomeSpaceMenu, helpMenu);
    final String os = System.getProperty ("os.name");
    if (os != null && os.startsWith ("Mac"))
      menuBar.useSystemMenuBarProperty ().set (true);

    VBox root = new VBox();
    Scene scene = new Scene(root, 1000, 600, Color.WHITE);
    root.getChildren().add(menuBar);

    String defaultGenome = "Human hg19";
    ObservableList<String> genomes = FXCollections.observableArrayList(defaultGenome, "chr1.fasta");
    ComboBox<String> genomeSelector = new ComboBox<String>(genomes);
    genomeSelector.setValue(defaultGenome);

    ObservableList<String> chromosomes = FXCollections.observableArrayList("chr1, chr2, chr3");
    ComboBox<String> chromosomeSelector = new ComboBox<String>(chromosomes);
    
    TextField jumpToTextField = new TextField();
    Label jumpToLabel = new Label("Go");
    jumpToLabel.setLabelFor(jumpToTextField);
    HBox jumpToPane = new HBox(jumpToTextField, jumpToLabel);
    jumpToPane.setAlignment(Pos.CENTER);

    Button homeButton = new Button("");
    homeButton.setId("homeButton");
    Button leftArrowButton = new Button("");
    leftArrowButton.setId("leftArrowButton");
    Button rightArrowButton = new Button("");
    rightArrowButton.setId("rightArrowButton");
    Button refreshScreenButton = new Button("");
    refreshScreenButton.setId("refreshScreenButton");
    Button regionToolButton = new Button("");
    regionToolButton.setId("regionToolButton");
    Button resizeToWindowButton = new Button("");
    resizeToWindowButton.setId("resizeToWindowButton");
    Button infoSelectButton = new Button("");
    infoSelectButton.setId("infoSelectButton");
    Button rulerButton = new Button("");
    rulerButton.setId("rulerButton");
    
    Slider zoomLevelSlider = new Slider();
    zoomLevelSlider.setShowTickMarks(true);
    zoomLevelSlider.setSnapToTicks(true);

    ToolBar toolbar = new ToolBar(genomeSelector, chromosomeSelector, jumpToPane,
        homeButton, leftArrowButton, rightArrowButton, refreshScreenButton, regionToolButton, resizeToWindowButton,
        infoSelectButton, rulerButton, zoomLevelSlider);
    root.getChildren().add(toolbar);

    ResizableCanvas resizableCanvas1 = new ResizableCanvas();
    Canvas canvas = resizableCanvas1.getCanvas();
    canvas.setHeight(500);
    canvas.setWidth(1000);
    
    ResizableCanvas resizableCanvas2 = new ResizableCanvas();
    Canvas canvas2 = resizableCanvas2.getCanvas();
    canvas2.setHeight(200);
    canvas2.setWidth(1000);
    
    SplitPane trackPane = new SplitPane(resizableCanvas1, resizableCanvas2);
    trackPane.setOrientation(Orientation.VERTICAL);
    trackPane.setDividerPositions(0.9);
    root.getChildren().add(trackPane);
    
    // Drawing experiment, nothing real
    GraphicsContext graphicsContext = canvas.getGraphicsContext2D();
    
    graphicsContext.setStroke(Color.RED);
    graphicsContext.setLineWidth(1);
    graphicsContext.setLineCap(StrokeLineCap.BUTT);
    graphicsContext.setLineDashes(10d, 5d, 15d, 20d);
    graphicsContext.setLineDashOffset(0);
    graphicsContext.strokeLine(10, 10, 200, 10);
    
    graphicsContext.setStroke(Color.GREEN);
    graphicsContext.setLineWidth(1);
    graphicsContext.setLineCap(StrokeLineCap.ROUND);
    graphicsContext.strokeLine(10, 30, 200, 30);
    
    graphicsContext.setStroke(Color.BLUE);
    graphicsContext.setLineWidth(1);
    graphicsContext.strokeLine(10, 50, 200, 50);

    GraphicsContext graphicsContext2 = canvas2.getGraphicsContext2D();
    
    graphicsContext2.setStroke(Color.PURPLE);
    graphicsContext2.setLineWidth(1);
    graphicsContext2.setLineCap(StrokeLineCap.BUTT);
    graphicsContext2.strokeLine(10, 10, 200, 10);

    scene.getStylesheets().add(getClass().getResource("experiment.css").toExternalForm());
    
    primaryStage.setScene(scene);
    primaryStage.show();
  }

  public static void main(String[] args) {
    launch(args);
  }
}
