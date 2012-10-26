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

import com.jidesoft.plaf.LookAndFeelFactory;
import jargs.gnu.CmdLineParser;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.event.GlobalKeyDispatcher;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.StringUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;


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

        initApplication();

        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ImageIcon icon = new ImageIcon(Main.class.getResource("mainframeicon.png"));
        if (icon != null) frame.setIconImage(icon.getImage());
        open(frame, args);

    }

    private static void initApplication() {
        DirectoryManager.initializeLog();
        log.info("Startup  " + Globals.applicationString());
        log.info("Java " + System.getProperty(Globals.JAVA_VERSION_STRING));
        log.info("Default User Directory: " + DirectoryManager.getUserDirectory());
        System.setProperty("http.agent", Globals.applicationString());

        Runtime.getRuntime().addShutdownHook(new ShutdownThread());

        updateTooltipSettings();

        // Anti alias settings.   TODO = Are these neccessary anymore ?
        System.setProperty("awt.useSystemAAFontSettings", "on");
        System.setProperty("swing.aatext", "true");


    }

    public static void updateTooltipSettings() {
        ToolTipManager.sharedInstance().setEnabled(true);
        final PreferenceManager prefMgr = PreferenceManager.getInstance();
        ToolTipManager.sharedInstance().setInitialDelay(prefMgr.getAsInt(PreferenceManager.TOOLTIP_INITIAL_DELAY));
        ToolTipManager.sharedInstance().setReshowDelay(prefMgr.getAsInt(PreferenceManager.TOOLTIP_RESHOW_DELAY));
        ToolTipManager.sharedInstance().setDismissDelay(prefMgr.getAsInt(PreferenceManager.TOOLTIP_DISMISS_DELAY));

    }


    /**
     * Open an IGV instance in the supplied Frame.
     *
     * @param frame
     */
    public static void open(Frame frame) {

        open(frame, new String[]{});
    }


    /**
     * Open an IGV instance in the supplied frame.
     *
     * @param frame
     * @param args  command-line arguments
     */
    public static void open(Frame frame, String[] args) {

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

        Main.IGVArgs igvArgs = new Main.IGVArgs(args);

        // Optional arguments
        if (igvArgs.getPropertyOverrides() != null) {
            PreferenceManager.getInstance().loadOverrides(igvArgs.getPropertyOverrides());
        }
        if (igvArgs.getDataServerURL() != null) {
            PreferenceManager.getInstance().overrideDataServerURL(igvArgs.getDataServerURL());
        }
        if (igvArgs.getGenomeServerURL() != null) {
            PreferenceManager.getInstance().overrideGenomeServerURL(igvArgs.getGenomeServerURL());
        }


        HttpUtils.getInstance().updateProxySettings();

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
        private String name = null;

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
            CmdLineParser.Option indexFileOption = parser.addStringOption('i', "indexFileURL");
            CmdLineParser.Option nameOption = parser.addStringOption('n', "name");

            try {
                parser.parse(args);
            } catch (CmdLineParser.IllegalOptionValueException e) {
                e.printStackTrace();  // This is not logged because the logger is not initialized yet.
            } catch (CmdLineParser.UnknownOptionException e) {
                e.printStackTrace();
            }
            propertyOverrides = (String) parser.getOptionValue(propertyFileOption);
            batchFile = (String) parser.getOptionValue(batchFileOption);
            port = (String) parser.getOptionValue(portOption);
            genomeId = (String) parser.getOptionValue(genomeOption);
            dataServerURL = (String) parser.getOptionValue(dataServerOption);
            genomeServerURL = (String) parser.getOptionValue(genomeServerOption);
            indexFile = (String) parser.getOptionValue(indexFileOption);
            name = (String) parser.getOptionValue(nameOption);

            String[] nonOptionArgs = parser.getRemainingArgs();
            if (nonOptionArgs != null && nonOptionArgs.length > 0) {
                String firstArg = StringUtils.decodeURL(nonOptionArgs[0]);  // TODO -- why is this url decoded?
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

// Alternative implementation
//                String firstArg = StringUtils.decodeURL(nonOptionArgs[0]);
//                firstArg=checkEqualsAndExtractParamter(firstArg);
//                if (firstArg != null && !firstArg.equals("ignore")) {
//                    log.info("Loading: " + firstArg);
//                    if (firstArg.endsWith(".xml") || firstArg.endsWith(".php") || firstArg.endsWith(".php3")
//                            || firstArg.endsWith(".session")) {
//                        sessionFile = firstArg;
//
//                    } else {
//                        dataFileString = firstArg;
//                    }
//                }
//
//                if (nonOptionArgs.length > 1) {
//                    // check if arg contains = for all args
//                    for (String arg: nonOptionArgs ) {
//                        arg = checkEqualsAndExtractParamter(arg);
//                        if (arg != null) locusString = arg;
//
//                    }
//
//                }
            }
        }

        private String checkEqualsAndExtractParamter(String arg) {
            if (arg == null) return null;
            int eq = arg.indexOf("=");
            if (eq > 0) {
                // we got a key=value
                String key = arg.substring(0, eq);
                String val = arg.substring(eq + 1);

                if (key.equalsIgnoreCase("server")) {
                    PreferenceManager.getInstance().put(PreferenceManager.IONTORRENT_SERVER, val);
                    log.info("Got server: " + key + "=" + val);
                    return null;
                } else if (key.equalsIgnoreCase("sessionURL") || key.equalsIgnoreCase("file")) {

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

        public String getName() {
            return name;
        }
    }

}
