package org.broad.igv.ui.javafx;

import htsjdk.samtools.seekablestream.SeekableStreamFactory;

import java.io.InputStream;
import java.util.List;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.DefaultExceptionHandler;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.Main;
import org.broad.igv.ui.ShutdownThread;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

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
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.Slider;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TextField;
import javafx.scene.control.ToolBar;
import javafx.scene.image.Image;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.StrokeLineCap;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

public class MainApplication extends Application {

  private static Logger log = Logger.getLogger(MainApplication.class);

  public MainApplication() {
  }

  @Override
  public void init() throws Exception {
    super.init();
    
    // TODO: refactor IGVArgs to be a stand-alone class or move it here as an inner class.
    // Note that JavaFX has a native Parameters class that we could use directly instead, though it's not an urgent
    // thing to switch.  We *do* need to use it go get the args array so we can get access at the instance level.
    // It shouldn't be a hard rewrite to switch and would remove the jargs dependency; it's just not important right now.
    List<String> args = getParameters().getRaw();
    final Main.IGVArgs igvArgs = new Main.IGVArgs(args);
    
    // Do this early
    if (igvArgs.igvDirectory != null) {
      // TODO: relocate the setIgvDirectory method to this class. 
        Main.setIgvDirectory(igvArgs);
    }
    
    // Application initialization.  Some of this could happen in the constructor, or in static blocks, etc,
    // but we'll need to re-evaluate where it should live (and whether its all still needed).

    // Win10 workaround?  Might not be needed since FileDialog will be a different widget.
    // Anti-aliasing workaround?  Also skipping for now..
    
    DirectoryManager.initializeLog();
    log.info("Startup  " + Globals.applicationString());
    log.info("Java " + System.getProperty(Globals.JAVA_VERSION_STRING));
    log.info("Default User Directory: " + DirectoryManager.getUserDirectory());
    log.info("OS: " + System.getProperty("os.name"));
    System.setProperty("http.agent", Globals.applicationString());

    // Add shutdown hook?  Can maybe happen in this class' stop() method (see discussion there).
    Runtime.getRuntime().addShutdownHook(new ShutdownThread());

    // - Tooltip management?
    // -- Apparently the delay setting won't be available until Java 9.  Discussed here, including a hack workaround:
    //    https://stackoverflow.com/questions/26854301/how-to-control-the-javafx-tooltips-delay
    // -- We'll skip this for now; can add it later.
    
    // TODO: relocate the checkVersion() method to this class
    Main.checkVersion();
    
    // Possible UI settings: in Swing it's necessary to wire up a lot of basic stuff (windowClose, etc).
    // Some or all of these may be unnecessary, or new ones might be required.  Also, these might need to
    // happen in the start() method, depending on what they do.
    // - windowClose handling
    // - enable tooltips
    // - L&F?  JavaFX has its own "look" / skin.  Need to figure out what we want here.
    // - HiDPI scaling.  This is important, but let's see if it's required.
    // - Keyboard dispatch

        // Optional arguments
    if (igvArgs.getPropertyOverrides() != null) {
        PreferencesManager.loadOverrides(igvArgs.getPropertyOverrides());
    }
    if (igvArgs.getDataServerURL() != null) {
        PreferencesManager.getPreferences().overrideDataServerURL(igvArgs.getDataServerURL());
    }
    if (igvArgs.getGenomeServerURL() != null) {
        PreferencesManager.getPreferences().overrideGenomeServerURL(igvArgs.getGenomeServerURL());
    }
    
    HttpUtils.getInstance().updateProxySettings();

    SeekableStreamFactory.setInstance(IGVSeekableStreamFactory.getInstance());

    // Load plug-in JARs.  Is this UI stuff?  Got to find out what these are about.  Suspect it's eps, but what else?
    RuntimeUtils.loadPluginJars();
  }

  private void checkLowMemory() {
    long mem = RuntimeUtils.getAvailableMemory();
    int MB = 1000000;
    
    if (mem < 400 * MB) {
        int mb = (int) (mem / MB);
        Alert alert = new Alert(AlertType.WARNING, "IGV is running with low available memory (" + mb + " mb)");
        alert.showAndWait();
    }
  }
  
  @Override
  public void start(Stage primaryStage) throws Exception {
    checkLowMemory();
    
    primaryStage.setTitle("JavaFX IGV");

    // This correctly loads the Frame (Application) icon, but it looks like this is not used by the Swing code
    InputStream iconStream = MainApplication.class.getResourceAsStream("/org/broad/igv/ui/mainframeicon.png");
    if (iconStream != null) {
      primaryStage.getIcons().add(new Image(iconStream));
    }

    // Here, launch a JavaFX version of the ui.IGV class.  All of the following needs to be moved into that class
    // or to components that it owns.  That is, once it launches then control should exit here (I think).
    // TODO: need the build a JavaFX version of this:
    // IGV.createInstance(stage).startUp(igvArgs);
    
    MenuItem loadFromFile = new MenuItem("Load from File ...");
    loadFromFile.setOnAction(new EventHandler<ActionEvent>() {
      
      @Override
      public void handle(ActionEvent actionEvent) {
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

    scene.getStylesheets().add(MainApplication.class.getResource("igv.css").toExternalForm());
    
    primaryStage.setScene(scene);
    primaryStage.show();
  }

  @Override
  public void stop() throws Exception {
    super.stop();
    
    // Application clean-up
    // Big question here is, how much can this take the place of ShutdownThread?
    // - I suspect that we still want that since it hooks into the JVM at a more basic level and doesn't rely on JavaFX
    // - Maybe we can handle basic stuff here, though
    
  }

  public static void main(String[] args) {
    Thread.setDefaultUncaughtExceptionHandler(new DefaultExceptionHandler());

    // Signal that we're running the JavaFX UI.
    Globals.IS_JAVAFX_UI = true;
    
    launch(args);
  }
}
