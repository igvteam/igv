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

import com.jidesoft.plaf.LookAndFeelFactory;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import jargs.gnu.CmdLineParser;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.Arrays;
import java.util.HashSet;

import static org.broad.igv.prefs.Constants.*;


/**
 * Utility class for launching IGV.  Provides a "main" method and an "open"  method for opening IGV in a supplied Frame.
 * <p/>
 * Note: The "open" methods must be executed on the event thread, for example
 * <p/>
 * public static void main(String[] args) {
 * EventQueue.invokeLater(new Runnable() {
 * public void run() {
 * Frame frame = new Frame();
 * org.broad.igv.ui.Main.open(frame);
 * }
 * );
 * }
 *
 * @author jrobinso
 * @date Feb 7, 2011
 */
public class Main {

    private static Logger log = Logger.getLogger(Main.class);

    /**
     * Launch an igv instance as a stand-alone application in its own Frame.
     *
     * @param args
     */
    public static void main(final String args[]) {

        Thread.setDefaultUncaughtExceptionHandler(new DefaultExceptionHandler());

        final Main.IGVArgs igvArgs = new Main.IGVArgs(args);

        // Do this early
        if (igvArgs.igvDirectory != null) {
            setIgvDirectory(igvArgs);
        }

        Runnable runnable = new Runnable() {
            public void run() {

                // This is a workaround for an internal JVM crash that was happening on Windows 10 (Creators Update).
                // TODO: remove when enough users have migrated to Java 8u141 or greater.
                // http://bugs.java.com/bugdatabase/view_bug.do?bug_id=JDK-8179014
                if (Globals.IS_WINDOWS && System.getProperty("os.name").contains("10")) {
                    UIManager.put("FileChooser.useSystemExtensionHiding", false);
                }

                String javaVersion = System.getProperty("java.version");
                if (javaVersion == null || !javaVersion.startsWith("1.8")) {
                    try {
                        System.out.println("Detected an unsupported Java version.  Java 8 is required by this release.");

                        if (!GraphicsEnvironment.isHeadless()) {
                            JOptionPane.showMessageDialog(null, "Detected an unsupported Java version.  Java 8 is required by this release.");
                        }
                    } finally {
                        System.exit(1);
                    }
                }

                initApplication();

                JFrame frame = new JFrame();
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                ImageIcon icon = new ImageIcon(Main.class.getResource("mainframeicon.png"));
                if (icon != null) frame.setIconImage(icon.getImage());
                open(frame, igvArgs);
            }
        };

        SwingUtilities.invokeLater(runnable);

    }

    private static void setIgvDirectory(IGVArgs igvArgs) {

        File dir = new File(igvArgs.igvDirectory);
        if (!dir.exists()) {

            // doesn't exist -- try to create it
            try {
                dir.mkdir();
            } catch (Exception e) {
                log.error("Error creating igv directory " + dir.getAbsolutePath(), e);
                return;
            }

            if (dir.isDirectory()) {
                if (dir.canWrite()) {
                    DirectoryManager.setIgvDirectory(dir);
                } else {
                    log.error("IGV directory '" + dir.getAbsolutePath() + "'is not writable");
                }
            } else {
                log.error("'" + dir.getAbsolutePath() + "' is not a directory");
            }
        } else {
            log.error("'" + dir.getAbsolutePath() + "' not found");
        }

    }

    private static void initApplication() {

        long mem = RuntimeUtils.getAvailableMemory();
        int MB = 1000000;
        if (mem < 400 * MB) {
            int mb = (int) (mem / MB);
            JOptionPane.showMessageDialog(null, "Warning: IGV is running with low available memory (" + mb + " mb)");
        }

        DirectoryManager.initializeLog();
        log.info("Startup  " + Globals.applicationString());
        log.info("Java " + System.getProperty(Globals.JAVA_VERSION_STRING));
        log.info("Default User Directory: " + DirectoryManager.getUserDirectory());
        log.info("OS: " + System.getProperty("os.name"));
        System.setProperty("http.agent", Globals.applicationString());

        Runtime.getRuntime().addShutdownHook(new ShutdownThread());

        updateTooltipSettings();

        // Anti alias settings.   TODO = Are these neccessary anymore ?
        System.setProperty("awt.useSystemAAFontSettings", "on");
        System.setProperty("swing.aatext", "true");

        checkVersion();


    }

    public static void updateTooltipSettings() {
        ToolTipManager.sharedInstance().setEnabled(true);
        final IGVPreferences prefMgr = PreferencesManager.getPreferences();
        ToolTipManager.sharedInstance().setInitialDelay(prefMgr.getAsInt(TOOLTIP_INITIAL_DELAY));
        ToolTipManager.sharedInstance().setReshowDelay(prefMgr.getAsInt(TOOLTIP_RESHOW_DELAY));
        ToolTipManager.sharedInstance().setDismissDelay(prefMgr.getAsInt(TOOLTIP_DISMISS_DELAY));

    }

    private static void checkVersion() {

        Runnable runnable = () -> {
            try {
                Version thisVersion = Version.getVersion(Globals.VERSION);
                if (thisVersion != null) {

                    final String serverVersionString = HttpUtils.getInstance().getContentsAsString(new URL(Globals.getVersionURL())).trim();
                    // See if user has specified to skip this update

                    final String skipString = PreferencesManager.getPreferences().get(SKIP_VERSION);
                    if(skipString != null) {
                        HashSet<String> skipVersion = new HashSet<>(Arrays.asList(skipString.split(",")));
                        if (skipVersion.contains(serverVersionString)) return;
                    }

                    Version serverVersion = Version.getVersion(serverVersionString.trim());
                    if (serverVersion == null) return;

                    if (thisVersion.lessThan(serverVersion)) {
                        log.info("A later version of IGV is available (" + serverVersionString + ")");
                    }
                } else if (Globals.VERSION.contains("3.0_beta") || Globals.VERSION.contains("snapshot")) {
                    HttpUtils.getInstance().getContentsAsString(new URL(Globals.getVersionURL())).trim();
                } else {
                    log.info("Unknown version: " + Globals.VERSION);
                }

            } catch (Exception e) {
                log.error("Error checking version", e);
            } finally {

            }
        };

        (new Thread(runnable)).start();

    }


    /**
     * Open an IGV instance in the supplied Frame.  This method used by unit tests
     *
     * @param frame
     */
    public static void open(Frame frame) {

        open(frame, new IGVArgs(new String[]{}));
    }


    /**
     * Open an IGV instance in the supplied frame.
     *
     * @param frame
     * @param igvArgs command-line arguments
     */
    public static void open(Frame frame, Main.IGVArgs igvArgs) {

        // Add a listener for the "close" icon, unless its a JFrame
        if (!(frame instanceof JFrame)) {
            frame.addWindowListener(new WindowAdapter() {
                @Override
                public void windowClosing(WindowEvent windowEvent) {
                    windowEvent.getComponent().setVisible(false);
                }
            });
        }

        // Turn on tooltip in case it was disabled for a temporary keyboard event, e.g. alt-tab
        frame.addWindowListener(new WindowAdapter() {

            @Override
            public void windowActivated(WindowEvent e) {
                ToolTipManager.sharedInstance().setEnabled(true);
            }

            @Override
            public void windowGainedFocus(WindowEvent windowEvent) {
                this.windowActivated(windowEvent);
            }
        });

        initializeLookAndFeel();

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
        IGV.createInstance(frame).startUp(igvArgs);

        // TODO Should this be done here?  Will this step on other key dispatchers?
        KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new GlobalKeyDispatcher());
    }

    private static void initializeLookAndFeel() {

        try {
            String lnf = UIManager.getSystemLookAndFeelClassName();
            UIManager.setLookAndFeel(lnf);

        } catch (Exception e) {
            e.printStackTrace();
        }

        double resolutionScale = Toolkit.getDefaultToolkit().getScreenResolution() / Globals.DESIGN_DPI;
        final IGVPreferences prefMgr = PreferencesManager.getPreferences();
        if (resolutionScale > 1.5) {
            if (prefMgr.getAsBoolean(SCALE_FONTS)) {
                FontManager.scaleFontSize(resolutionScale);
            } else if (prefMgr.hasExplicitValue(DEFAULT_FONT_SIZE)) {
                int fs = prefMgr.getAsInt(DEFAULT_FONT_SIZE);
                FontManager.updateSystemFontSize(fs);
            }
        }


        if (Globals.IS_LINUX) {
            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
                UIManager.put("JideSplitPane.dividerSize", 5);
                UIManager.put("JideSplitPaneDivider.background", Color.darkGray);

            } catch (Exception exception) {
                exception.printStackTrace();
            }

        }
        // Todo -- what does this do?
        LookAndFeelFactory.installJideExtension();
    }


    /**
     * Class to represent the IGV version. The version numbering follows Microsoft conventions.  This class is public
     * to allow unit tests.
     * <p/>
     * Example internal version string:  2.3.27 (31)02/18/2014 11:42 PM
     * Example server version string:  2.3.27
     */
    public static class Version {

        private int major;
        private int minor;
        private int build;

        public static Version getVersion(String versionString) {
            String[] tokens = versionString.split("\\.");
            if (tokens.length < 2) {
                return null;   // Unknown version
            } else {
                try {
                    int major = Integer.parseInt(tokens[0]);
                    int minor = Integer.parseInt(tokens[1]);
                    int build = tokens.length <= 2 ? 0 : Integer.parseInt(tokens[2]);
                    return new Version(major, minor, build);

                } catch (NumberFormatException e) {
                    log.error("Error parsing version string: " + versionString);
                    return null;
                }
            }
        }

        private Version(int major, int minor, int build) {
            this.major = major;
            this.minor = minor;
            this.build = build;
        }

        public int getMajor() {
            return major;
        }

        public int getMinor() {
            return minor;
        }

        public int getBuild() {
            return build;
        }

        public boolean lessThan(Version anotherVersion) {
            if (anotherVersion.major > this.major) return true;
            else if (anotherVersion.major < this.major) return false;
            else if (anotherVersion.minor > this.minor) return true;
            else if (anotherVersion.minor < this.minor) return false;
            else return anotherVersion.build > this.build;
        }
    }

    /**
     * Class to encapsulate IGV command line arguments.
     */
    static public class IGVArgs {
        private String batchFile = null;
        private String sessionFile = null;
        private String dataFileString = null;
        private String locusString = null;
        private String propertyOverrides = null;
        private String genomeId = null;
        private String port = null;
        private String dataServerURL = null;
        private String genomeServerURL = null;
        private String indexFile = null;
        private String coverageFile = null;
        private String name = null;
        public String igvDirectory = null;
        public String forceVersion = null;

        IGVArgs(String[] args) {
            if (args != null) {
                parseArgs(args);
            }
        }

        /**
         * Parse arguments.  All arguments are optional,  a full set of arguments are
         * firstArg  locusString  -b batchFile -p preferences
         */
        private void parseArgs(String[] args) {
            CmdLineParser parser = new CmdLineParser();
            CmdLineParser.Option propertyFileOption = parser.addStringOption('o', "preferences");
            CmdLineParser.Option batchFileOption = parser.addStringOption('b', "batch");
            CmdLineParser.Option portOption = parser.addStringOption('p', "port");
            CmdLineParser.Option genomeOption = parser.addStringOption('g', "genome");
            CmdLineParser.Option dataServerOption = parser.addStringOption('d', "dataServerURL");
            CmdLineParser.Option genomeServerOption = parser.addStringOption('u', "genomeServerURL");
            CmdLineParser.Option indexFileOption = parser.addStringOption('i', "indexFile");
            CmdLineParser.Option coverageFileOption = parser.addStringOption('c', "coverageFile");
            CmdLineParser.Option nameOption = parser.addStringOption('n', "name");
            CmdLineParser.Option igvDirectoryOption = parser.addStringOption("igvDirectory");
            CmdLineParser.Option forceVersionOption = parser.addStringOption("forceVersion");

            try {
                parser.parse(args);
            } catch (CmdLineParser.IllegalOptionValueException e) {
                e.printStackTrace();  // This is not logged because the logger is not initialized yet.
            } catch (CmdLineParser.UnknownOptionException e) {
                e.printStackTrace();
            }
            propertyOverrides = getDecodedValue(parser, propertyFileOption);
            batchFile = getDecodedValue(parser, batchFileOption);
            port = (String) parser.getOptionValue(portOption);
            genomeId = (String) parser.getOptionValue(genomeOption);
            dataServerURL = getDecodedValue(parser, dataServerOption);
            genomeServerURL = getDecodedValue(parser, genomeServerOption);
            name = (String) parser.getOptionValue(nameOption);

            String indexFilePath = (String) parser.getOptionValue(indexFileOption);
            if (indexFilePath != null) {
                indexFile = maybeDecodePath(indexFilePath);
            }

            String coverageFilePath = (String) parser.getOptionValue(coverageFileOption);
            if (coverageFilePath != null) {
                coverageFile = maybeDecodePath(coverageFilePath);
            }


            String igvDirectoryPath = (String) parser.getOptionValue(igvDirectoryOption);
            if (igvDirectoryPath != null) {
                igvDirectory = maybeDecodePath(igvDirectoryPath);
            }

            String forceVersion = (String) parser.getOptionValue(forceVersionOption);
            if(forceVersion != null) {
                Globals.VERSION = forceVersion;
            }

            String[] nonOptionArgs = parser.getRemainingArgs();

            if (nonOptionArgs != null && nonOptionArgs.length > 0) {

                String firstArg = maybeDecodePath(nonOptionArgs[0]);

                if (firstArg != null && !firstArg.equals("ignore")) {
                    log.info("Loading: " + firstArg);
                    if (firstArg.endsWith(".xml") || firstArg.endsWith(".php") || firstArg.endsWith(".php3")
                            || firstArg.endsWith(".session")) {
                        sessionFile = firstArg;
                    } else {
                        dataFileString = firstArg;
                    }
                }
                if (nonOptionArgs.length > 1) {
                    locusString = nonOptionArgs[1];
                }


            }
        }

        private String maybeDecodePath(String path) {

            if (FileUtils.isRemote(path)) {
                return URLDecoder.decode(path);
            } else {
                if (FileUtils.resourceExists(path)) {
                    return path;
                } else {
                    return URLDecoder.decode(path);
                }
            }
        }

        private String getDecodedValue(CmdLineParser parser, CmdLineParser.Option option) {

            String value = (String) parser.getOptionValue(option);
            if (value == null) return null;

            try {
                return URLDecoder.decode(value, "UTF-8");
            } catch (UnsupportedEncodingException e) {
                log.error(e);
                return value;
            }
        }

        private String checkEqualsAndExtractParamter(String arg) {
            if (arg == null) return null;
            int eq = arg.indexOf("=");
            if (eq > 0) {
                // we got a key=value
                String key = arg.substring(0, eq);
                String val = arg.substring(eq + 1);

                if (key.equalsIgnoreCase("sessionURL") || key.equalsIgnoreCase("file")) {

                    if (val.endsWith(".xml") || val.endsWith(".php") || val.endsWith(".php3")
                            || val.endsWith(".session")) {
                        log.info("Got session: " + key + "=" + val);
                        sessionFile = val;

                    } else {
                        log.info("Got dataFileString: " + key + "=" + val);
                        dataFileString = val;
                    }
                    return null;
                } else if (key.equalsIgnoreCase("locus") || key.equalsIgnoreCase("position")) {
                    log.info("Got locus: " + key + "=" + val);
                    locusString = val;
                    return null;
                } else {
                    log.info("Currently not handled: " + key + "=" + val);
                    return null;
                }

            }
            return arg;
        }

        public String getBatchFile() {
            return batchFile;
        }

        public String getSessionFile() {
            return sessionFile;
        }

        public String getDataFileString() {
            return dataFileString;
        }

        public String getLocusString() {
            return locusString;
        }

        public String getPropertyOverrides() {
            return propertyOverrides;
        }

        public String getGenomeId() {
            return genomeId;
        }

        public String getPort() {
            return port;
        }

        public String getDataServerURL() {
            return dataServerURL;
        }

        public String getGenomeServerURL() {
            return genomeServerURL;
        }

        public String getIndexFile() {
            return indexFile;
        }

        public String getCoverageFile() {
            return coverageFile;
        }

        public String getName() {
            return name;
        }

    }

}
