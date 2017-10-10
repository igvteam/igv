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

import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import javafx.application.Application;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.stage.Stage;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.Session;
import org.broad.igv.ui.DefaultExceptionHandler;
import org.broad.igv.ui.Main;
import org.broad.igv.ui.ShutdownThread;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.util.List;

// Intended as the rough equivalent of the Main class of the Swing UI.  Work in progress.
public class MainApplication extends Application {

    private static Logger log = Logger.getLogger(MainApplication.class);

    public MainApplication() {
    }

    private Main.IGVArgs igvArgs;

    @Override
    public void init() throws Exception {
        super.init();

        // TODO: refactor IGVArgs to be a stand-alone class or move it here as an inner
        // class.
        // Note that JavaFX has a native Parameters class that we could use directly
        // instead, though it's not an urgent
        // thing to switch. We *do* need to use it go get the args array so we can get
        // access at the instance level.
        // It shouldn't be a hard rewrite to switch and would remove the jargs
        // dependency; it's just not important right now.
        List<String> args = getParameters().getRaw();
        igvArgs = new Main.IGVArgs(args);

        // Do this early
        if (igvArgs.igvDirectory != null) {
            // TODO: relocate the setIgvDirectory method to this class.
            Main.setIgvDirectory(igvArgs);
        }

        // Application initialization. Some of this could happen in the constructor, or
        // in static blocks, etc,
        // but we'll need to re-evaluate where it should live (and whether its all still
        // needed).

        // Win10 workaround? Might not be needed since FileDialog will be a different
        // widget.
        // Anti-aliasing workaround? Also skipping for now..

        DirectoryManager.initializeLog();
        log.info("Startup  " + Globals.applicationString());
        log.info("Java " + System.getProperty(Globals.JAVA_VERSION_STRING));
        log.info("Default User Directory: " + DirectoryManager.getUserDirectory());
        log.info("OS: " + System.getProperty("os.name"));
        System.setProperty("http.agent", Globals.applicationString());

        // Add shutdown hook? Can maybe happen in this class' stop() method (see
        // discussion there).
        Runtime.getRuntime().addShutdownHook(new ShutdownThread());

        // - Tooltip management?
        // -- Apparently the delay setting won't be available until Java 9. Discussed
        // here, including a hack workaround:
        // https://stackoverflow.com/questions/26854301/how-to-control-the-javafx-tooltips-delay
        // -- We'll skip this for now; can add it later.

        // TODO: relocate the checkVersion() method to this class
        Main.checkVersion();

        // Possible UI settings: in Swing it's necessary to wire up a lot of basic stuff
        // (windowClose, etc).
        // Some or all of these may be unnecessary, or new ones might be required. Also,
        // these might need to
        // happen in the start() method, depending on what they do.
        // - windowClose handling
        // - enable tooltips
        // - L&F? JavaFX has its own "look" / skin. Need to figure out what we want
        // here.
        // - HiDPI scaling. This is important, but let's see if it's required.
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

        // TODO: port these to JavaFX (Swing refs)
        // These may more properly belong elsewhere...
        GenomeManager genomeManager = GenomeManager.getInstance();
        Session session = new Session(null);

        // Here, launch a JavaFX version of the ui.IGV class.
        // TODO: working on a JavaFX version of this
        IGVStageBuilder.buildStage(primaryStage);

        primaryStage.show();
    }

    @Override
    public void stop() throws Exception {
        super.stop();

        // Application clean-up
        // Big question here is, how much can this take the place of ShutdownThread?
        // - I suspect that we still want that since it hooks into the JVM at a more
        // basic level and doesn't rely on JavaFX
        // - Maybe we can handle basic stuff here, though

    }

    public static void main(String[] args) {
        Thread.setDefaultUncaughtExceptionHandler(new DefaultExceptionHandler());

        // Signal that we're running the JavaFX UI.
        Globals.IS_JAVAFX_UI = true;

        launch(args);
    }
}
