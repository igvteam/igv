/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.ui;

import com.jidesoft.plaf.LookAndFeelFactory;
import jargs.gnu.CmdLineParser;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.event.GlobalKeyDispatcher;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.IGVHttpUtils;

import javax.swing.*;
import java.awt.*;


/**
 * @author jrobinso
 * @date Feb 7, 2011
 */
public class Main {

    private static Logger log = Logger.getLogger(Main.class);

 
    public static void main(final String args[]) {


        //RepaintManager.setCurrentManager(new TracingRepaintManager());
        log.info("Startup");
        log.info(Globals.applicationString());

        System.setProperty("http.agent", Globals.applicationString());

        FileUtils.addRollingAppenderToRootLogger();
        log.info("Default User Directory: " + Globals.getUserDirectory());

        Thread.setDefaultUncaughtExceptionHandler(new DefaultExceptionHandler());


        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                com.jidesoft.utils.Lm.verifyLicense("The Broad Institute, MIT", "Gene Pattern",
                        "D.DQSR7z9m6fxL1IqWZ6svQFmE6vj3Q");

                // Set look and feel
                if (!Globals.IS_MAC) {
                    try {
                        for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
                            if ("Nimbus".equals(info.getName())) {
                                UIManager.setLookAndFeel(info.getClassName());
                                break;
                            }
                        }
                    }
                    catch (Exception e) {
                        log.error("Error installing look and feel", e);
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

                } // Todo -- what does this do?
                LookAndFeelFactory.installJideExtension();

                IGVMainFrame frame = null;
                JWindow splashScreen = null;
                try {

                    frame = IGVMainFrame.getInstance();

                    IGVHttpUtils.updateProxySettings();

                    frame.startUp(args);

                    KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new GlobalKeyDispatcher());

                } catch (Exception e) {

                    log.error("Fatal application error!", e);
                    System.exit(-1);
                } finally {
                    if (splashScreen != null) {
                        splashScreen.setVisible(false);
                    }

                }
            }
        });
    }



    /**
     * Class to encapsulate IGV command line arguments.
     */
    static class IGVArgs {
        private String batchFile = null;
        private String sessionFile = null;
        private String dataFileString = null;
        private String locusString = null;
        private String propertyFile = null;
        private String genomeId = null;
        private String port = null;
        private String dataServerURL = null;
        private String genomeServerURL = null;

        IGVArgs(String[] args) {
            parseArgs(args);
        }

        /**
         * Parse arguments.  All arguments are optional,  a full set of arguments are
         * firstArg  locusString  -b batchFile -p preferences
         */
        private void parseArgs(String[] args) {
            CmdLineParser parser = new CmdLineParser();
            CmdLineParser.Option propertyFileOption = parser.addStringOption('p', "preferences");
            CmdLineParser.Option batchFileOption = parser.addStringOption('b', "batch");
            CmdLineParser.Option portOption = parser.addStringOption('p', "port");
            CmdLineParser.Option genomeOption = parser.addStringOption('g', "genome");
            CmdLineParser.Option dataServerOption = parser.addStringOption('d', "dataServerURL");
            CmdLineParser.Option genomeServerOption = parser.addStringOption('u', "genomeServerURL");

            try {
                parser.parse(args);
            } catch (CmdLineParser.IllegalOptionValueException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (CmdLineParser.UnknownOptionException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            propertyFile = (String) parser.getOptionValue(propertyFileOption);
            batchFile = (String) parser.getOptionValue(batchFileOption);
            port = (String) parser.getOptionValue(portOption);
            genomeId = (String) parser.getOptionValue(genomeOption);
            dataServerURL = (String) parser.getOptionValue(dataServerOption);
            genomeServerURL = (String) parser.getOptionValue(genomeServerOption);

            String[] nonOptionArgs = parser.getRemainingArgs();
            if (nonOptionArgs != null && nonOptionArgs.length > 0) {
                String firstArg = nonOptionArgs[0];
                if (!firstArg.equals("ignore")) {
                    if (firstArg.endsWith("xml")) {
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

        public String getPropertyFile() {
            return propertyFile;
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
    }
}
