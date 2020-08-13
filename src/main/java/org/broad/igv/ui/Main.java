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
import com.sanityinc.jargs.CmdLineParser;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.google.OAuthUtils;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.broad.igv.util.stream.IGVUrlHelperFactory;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileWriter;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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

        htsjdk.tribble.util.ParsingUtils.setURLHelperFactory(IGVUrlHelperFactory.getInstance());

        try {
            OAuthUtils.getInstance();  // Initialize oauth
        } catch (Exception e) {
            log.error("Warning: Error fetching oAuth properties: " + e.getMessage());
        }

        final Main.IGVArgs igvArgs = new Main.IGVArgs(args);

        // Do this early
        if (igvArgs.igvDirectory != null) {
            setIgvDirectory(igvArgs);
        }
        checkDotIgvDirectory();

        Runnable runnable = () -> {
            if (Globals.IS_WINDOWS && System.getProperty("os.name").contains("10")) {
                UIManager.put("FileChooser.useSystemExtensionHiding", false);
            }

            DesktopIntegration.verifyJavaPlatform();
            initApplication();

            JFrame frame = new JFrame();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            ImageIcon icon = new ImageIcon(Main.class.getResource("mainframeicon.png"));
            if (icon != null) frame.setIconImage(icon.getImage());
            open(frame, igvArgs);
        };

        SwingUtilities.invokeLater(runnable);

    }

    private static void setIgvDirectory(IGVArgs igvArgs) {

        File dir = new File(igvArgs.igvDirectory);
        if (!dir.exists()) {

            // doesn't exist -- try to create it
            try {
                if (dir.mkdir()) {
                    DirectoryManager.setIgvDirectory(dir);
                } else {
                    log.error("Unable to create igv directory " + dir.getAbsolutePath());
                }
            } catch (Exception e) {
                log.error("Error creating igv directory " + dir.getAbsolutePath(), e);
                return;
            }
        } else if (dir.isDirectory()) {
            if (dir.canWrite()) {
                DirectoryManager.setIgvDirectory(dir);
            } else {
                log.error("IGV directory '" + dir.getAbsolutePath() + "'is not writable");
            }
        } else {
            log.error("'" + dir.getAbsolutePath() + "' is not a directory");
        }
    }

    private static void checkDotIgvDirectory() {
        // Check if the .igv directory exists and create it if not.  This is a config
        // file with a known name and location, not intended to be moved by the user.
        // At present, this is only used by the launcher scripts and not the Java code.
        String userHome = System.getProperty("user.home");
        File dir = new File(userHome, ".igv");
        if (!dir.exists()) {
            // doesn't exist -- try to create it
            try {
                dir.mkdir();
            } catch (Exception e) {
                // Ignore the mkdir failure.  It's not necessary to even report this.
                // We'll proceed without it.
                return;
            }
        }
        
        // Also check if the java_arguments file exists and create it if not.  This is
        // likewise only used by the launcher scripts.  We create it here as a user
        // convenience.  Note that we skip it if ~/.igv is not a directory.
        if (dir.isDirectory()) {
            File argsFile = new File(dir, "java_arguments");
            if (!argsFile.exists()) {
                // doesn't exist -- try to create it
                try {
                    FileWriter argsFileWriter = new FileWriter(argsFile);
                    try {
                        argsFileWriter.append("# See https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt for tips on using this file.");
                        argsFileWriter.append(System.lineSeparator());
                        argsFileWriter.append("# Uncomment the following line for an 8 GB memory spec for IGV.");
                        argsFileWriter.append(System.lineSeparator());
                        argsFileWriter.append("# -Xmx8G");
                        argsFileWriter.append(System.lineSeparator());
                    } finally {
                        argsFileWriter.close();
                    }
                } catch (Exception e) {
                    // As above, ignore the write failure.
                }
            }
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
        log.info("Java " + System.getProperty(Globals.JAVA_VERSION_STRING)
                + " (build " + System.getProperty("java.vm.version")
                + ") " + System.getProperty("java.version.date", ""));
        log.info("Java Vendor: " + System.getProperty("java.vendor")
                + " " + System.getProperty("java.vendor.url", ""));
        log.info("JVM: " + System.getProperty("java.vm.name", "")
                + " " + System.getProperty("java.vendor.version", "")
                + "   " + System.getProperty("java.compiler", ""));
        log.info("Default User Directory: " + DirectoryManager.getUserDirectory());
        log.info("OS: " + System.getProperty("os.name") + " " + System.getProperty("os.version")
                + " " + System.getProperty("os.arch"));
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
                    if (skipString != null) {
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
                // ignore
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
        private java.util.List<String> dataFileStrings = null;
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
        public Collection<String> httpHeader = null;

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
            CmdLineParser.Option locusOption = parser.addStringOption('l', "locus");
            CmdLineParser.Option igvDirectoryOption = parser.addStringOption("igvDirectory");
            CmdLineParser.Option versionOption = parser.addBooleanOption("version");
            CmdLineParser.Option helpOption = parser.addBooleanOption("help");
            CmdLineParser.Option headerOption = parser.addStringOption('H', "header");


            try {
                parser.parse(args);
            } catch (Exception e) {
                e.printStackTrace();  // This is not logged because the logger is not initialized yet.
                JOptionPane.showMessageDialog(null, "Error parsing command line arguments: " + e.getMessage());
                return;
            }

            propertyOverrides = getDecodedValue(parser, propertyFileOption);
            batchFile = getDecodedValue(parser, batchFileOption);
            port = (String) parser.getOptionValue(portOption);
            genomeId = (String) parser.getOptionValue(genomeOption);
            dataServerURL = getDecodedValue(parser, dataServerOption);
            genomeServerURL = getDecodedValue(parser, genomeServerOption);
            name = (String) parser.getOptionValue(nameOption);
            locusString = (String) parser.getOptionValue(locusOption);
            httpHeader = parser.getOptionValues(headerOption);

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

            String[] nonOptionArgs = parser.getRemainingArgs();

            if (parser.getOptionValue(versionOption) != null) {
                System.out.println(Globals.VERSION);
                System.exit(0);
            }

            if(parser.getOptionValue(helpOption) != null) {
                printHelp();
                System.exit(0);
            }

            // The Mac app launcher sometimes inserts "" into the command line.  Filter empty strings
            if (nonOptionArgs.length > 0) {
                nonOptionArgs = Arrays.stream(nonOptionArgs).filter(s -> !s.isEmpty()).toArray(String[]::new);
            }

            if (nonOptionArgs != null && nonOptionArgs.length > 0) {

                if (locusString != null || nonOptionArgs.length > 2) {

                    decodeNonOptionArgs(nonOptionArgs);

                } else if (nonOptionArgs.length == 1) {

                    decodeNonOptionArgsLegacy(nonOptionArgs);

                } else {  //nonOptionArgs.length == 2

                    // If the second argument is not a file or URL we assume its a locus => legacy
                    String secondArg = maybeDecodePath(nonOptionArgs[1]);
                    if (FileUtils.isRemote(secondArg) || (new File(secondArg)).exists()) {
                        decodeNonOptionArgs(nonOptionArgs);
                    } else {
                        decodeNonOptionArgsLegacy(nonOptionArgs);
                    }
                }
            }
        }

        private void decodeNonOptionArgsLegacy(String[] nonOptionArgs) {

            dataFileStrings = new ArrayList<>();

            String firstArg = maybeDecodePath(nonOptionArgs[0]);

            if (firstArg != null) {
                log.info("Loading: " + firstArg);
                if (firstArg.endsWith(".xml") || firstArg.endsWith(".php") || firstArg.endsWith(".php3")
                        || firstArg.endsWith(".session")) {
                    sessionFile = firstArg;
                } else {
                    String[] paths = firstArg.split(",");
                    for (String p : paths) {
                        dataFileStrings.add(p);
                    }
                }
            }
            if (nonOptionArgs.length == 2) {

                locusString = nonOptionArgs[1];
            }
        }

        /**
         * IGV version > 2.5 => non-optional arguments are files
         *
         * @param nonOptionArgs
         */
        private void decodeNonOptionArgs(String[] nonOptionArgs) {
            dataFileStrings = new ArrayList<>();
            for (String arg : nonOptionArgs) {
                dataFileStrings.add(maybeDecodePath(arg));
            }

        }

        private String maybeDecodePath(String path) {

            if ((new File(path)).exists()) {
                return path;
            } else if (FileUtils.isRemote(path)) {
                return URLDecoder.decode(path);
            } else if (path.contains("%2C") || path.contains("%3F") || path.contains("%2B") || path.contains("%2F")) {
                return URLDecoder.decode(path);
            } else {
                return path;
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

        public String getBatchFile() {
            return batchFile;
        }

        public String getSessionFile() {
            return sessionFile;
        }

        public java.util.List<String> getDataFileStrings() {
            return dataFileStrings;
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

        public Collection<String> getHttpHeader() {
            return httpHeader;
        }

        private void printHelp() {
            System.out.println("Command line options:");
            System.out.println("Space delimited list of data files to load");
            System.out.println("--preferences, -o  Path or URL to a preference property file");
            System.out.println("--batch. -b  Path or url to a batch command file");
            System.out.println("--port, -p  IGV command port number (defaults to 60151)");
            System.out.println("--genome, -g  Genome ID (e.g hg19) or path or url to .genome or indexed fasta file");
            System.out.println("--dataServerURL, -d  Path or url to a data server registry file");
            System.out.println("--genomeServerURL, -u  Path or url to a genome server registry file");
            System.out.println("--indexFile, -i  Index file or comma delimited list of index files corresponding to data files");
            System.out.println("--coverageFile, -c  Coverage file or comma delimited list of coverage files corresponding to data files");
            System.out.println("--name, -n  Name or comma-delimited list of names for tracks corresponding to data files");
            System.out.println("--locus, -l  Initial locus");
            System.out.println("--header, -H http header to include with all requests for list of data files");
            System.out.println("--igvDirectory Path to the local igv directory.  Defaults to <user home>/igv");
            System.out.println("--version  Print the IGV version and exit");
            System.out.println("--help Print this output and exit");
        }
    }
}
