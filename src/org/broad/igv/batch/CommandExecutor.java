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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.batch;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.TrackManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.SnapshotUtilities;
import org.broad.igv.util.*;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class CommandExecutor {

    private static Logger log = Logger.getLogger(CommandExecutor.class);

    private File snapshotDirectory;


    private List<String> getArgs(String[] tokens) {
        List<String> args = new ArrayList(tokens.length);
        for (String s : tokens) {
            if (s.trim().length() > 0) {
                args.add(s.trim());
            }
        }
        return args;
    }

    public String execute(String command) {

        List<String> args = getArgs(StringUtils.breakQuotedString(command, ' ').toArray(new String[]{}));

        String result = "OK";
        final IGV mainFrame = IGV.getFirstInstance();

        System.out.println();
        log.debug("Executing: " + command);
        try {
            if (args.size() > 0) {

                String cmd = args.get(0).toLowerCase();
                String param1 = args.size() > 1 ? args.get(1) : null;
                String param2 = args.size() > 2 ? args.get(2) : null;
                String param3 = args.size() > 3 ? args.get(3) : null;

                if (cmd.equals("echo")) {
                    result = cmd;
                } else if (cmd.equals("goto")) {
                    result = goto1(param1);
                } else if (cmd.equals("snapshotdirectory")) {
                    result = setSnapshotDirectory(param1);

                } else if (cmd.equals("snapshot")) {
                    String filename = param1;
                    createSnapshot(filename);

                } else if ((cmd.equals("loadfile") || cmd.equals("load")) && param1 != null) {
                    result = load(param1);
                } else if (cmd.equals("hget") && args.size() > 3) {
                    result = hget(param1, param2, param3);
                } else if (cmd.equals("genome") && args.size() > 1) {
                    result = genome(param1);
                } else if (cmd.equals("new") || cmd.equals("reset") || cmd.equals("clear")) {
                    mainFrame.createNewSession(null);
                } else if (cmd.equals("sort")) {
                    sort(param1, param2, param3);
                } else if (cmd.equals("collapse")) {
                    String trackName = param1 == null ? null : param1.replace("\"", "").replace("'", "");
                    collapse(trackName);
                } else if (cmd.equals("expand")) {
                    String trackName = param1 == null ? null : param1.replace("\"", "").replace("'", "");
                    expand(trackName);
                } else if (cmd.equals("tweakdivider")) {
                    IGV.getFirstInstance().tweakPanelDivider();
                } else if (cmd.equals("maxpanelheight") && param1 != null) {
                    return setMaxPanelHeight(param1);
                } else if (cmd.equals("exit")) {
                    System.exit(0);
                } else {
                    log.error("UNKOWN COMMAND: " + command);
                    return "UNKOWN COMMAND: " + command;
                }
            } else {
                return "Empty command string";
            }
            IGV.getFirstInstance().doRefresh();

            if (RuntimeUtils.getAvailableMemoryFraction() < 0.5) {
                log.debug("Clearing caches");
                LRUCache.clearCaches();
            }
            log.debug("Finished execution: " + command + "  sleeping ....");
            Thread.sleep(2000);
            log.debug("Finished sleeping");

        } catch (Exception e) {
            log.error("Could not Parse Command", e);
            return "ERROR Could not Parse Command: " + e.toString();
        }
        log.info(result);

        return result;
    }

    private String setMaxPanelHeight(String param1) {
        try {
            Integer h = Integer.parseInt(param1.trim());
            SnapshotUtilities.MAX_PANEL_HEIGHT = h;
            return "OK";
        }
        catch (NumberFormatException e) {
            return "ERROR - max panel height value ('" + param1 + ".) must be a number";
        }
    }

    private String genome(String param1) {
        if (param1 == null) {
            return "ERROR missing genome parameter";
        }
        String result;
        String genomeID = param1;
        IGV.getFirstInstance().selectGenomeFromList(genomeID);
        result = "OK";
        return result;
    }

    private String hget(String param1, String param2, String param3) throws IOException {
        String result;
        String fileString = param1;
        String locusString = param2;
        String mergeValue = param3;
        boolean merge = mergeValue != null && mergeValue.equalsIgnoreCase("true");
        result = loadFiles(fileString, locusString, merge);
        return result;
    }

    private String load(String param1) throws IOException {
        if (param1 == null) {
            return "ERROR: missing path parameter";
        }
        String fileString = param1.replace("\"", "").replace("'", "");
        return loadFiles(fileString, null, true);
    }

    private String setSnapshotDirectory(String param1) {
        if (param1 == null) {
            return "ERROR: missing directory parameter";
        }

        String result;
        File parentDir = new File(param1);
        if (parentDir.exists()) {
            snapshotDirectory = parentDir;
            result = "OK";
        } else {
            parentDir.mkdir();
            if (parentDir.exists()) {
                snapshotDirectory = parentDir;
                result = "OK";
            } else {

                result = "ERROR: directory: " + param1 + " does not exist";
            }
        }
        return result;
    }

    private String goto1(String param1) {
        if (param1 == null) {
            return "ERROR: missing locus parameter";
        }
        String locus = param1;
        IGV.getFirstInstance().goToLocus(locus);
        return "OK";
    }

    private void collapse(String trackName) {
        if (trackName == null) {
            IGV.getFirstInstance().getTrackManager().collapseTracks();
        } else {
            IGV.getFirstInstance().getTrackManager().collapseTrack(trackName);
        }
        IGV.getFirstInstance().repaintDataPanels();
    }


    private void expand(String trackName) {
        if (trackName == null) {
            IGV.getFirstInstance().getTrackManager().expandTracks();
        } else {
            IGV.getFirstInstance().getTrackManager().expandTrack(trackName);
        }
        IGV.getFirstInstance().repaintDataPanels();
    }


    private void sort(String sortArg, String locusString, String param3) {
        TrackManager tm = IGV.getFirstInstance().getTrackManager();
        RegionScoreType regionSortOption = getRegionSortOption(sortArg);
        if (regionSortOption != null) {
            RegionOfInterest roi = null;
            if (locusString != null) {
                Locus locus = new Locus(locusString);
                if (locus.isValid()) {
                    roi = new RegionOfInterest(locus.getChr(), locus.getStart(), locus.getEnd(), "");
                }
            }
            tm.sortByRegionScore(roi, regionSortOption, FrameManager.getDefaultFrame());

        } else {
            Double location = null;
            if (param3 != null) {
                try {
                    location = new Double(param3.replace(",", ""));
                }
                catch (NumberFormatException e) {
                    log.info("Unexpected sort location argument (expected number): " + param3);
                }
            }
            if (location == null) {
                tm.sortAlignmentTracks(getAlignmentSortOption(sortArg));
            } else {
                tm.sortAlignmentTracks(getAlignmentSortOption(sortArg), location);
            }

        }
        IGV.getFirstInstance().repaintDataPanels();
    }

    private String loadFiles(final String fileString, final String locus, final boolean merge) throws IOException {

        log.debug("Run load files");
        WaitCursorManager.CursorToken token = null;
        try {
            token = WaitCursorManager.showWaitCursor();

            String[] files = fileString.split(",");
            List<ResourceLocator> fileLocators = new ArrayList<ResourceLocator>();
            List<String> sessionPaths = new ArrayList<String>();

            if (!merge) {
                IGV.getFirstInstance().createNewSession(null);
            }

            for (String f : files) {
                if (f.endsWith(".xml")) {
                    sessionPaths.add(f);
                } else {
                    ResourceLocator rl = new ResourceLocator(f);
                    fileLocators.add(rl);
                }
            }

            for (String sessionPath : sessionPaths) {
                InputStream is = null;
                try {
                    is = ParsingUtils.openInputStream(new ResourceLocator(sessionPath));
                    IGV.getFirstInstance().doRestoreSession(is, sessionPath, locus, merge);
                } finally {
                    if (is != null) {
                        try {
                            is.close();
                        } catch (IOException e) {
                            log.error(e.getMessage(), e);
                        }
                    }
                }
            }

            //TODO Find a better way to get the message back up to the socket.
            IGV.getFirstInstance().loadTracks(fileLocators);

            if (locus != null) {
                IGV.getFirstInstance().goToLocus(locus);
            }
        }
        finally {
            WaitCursorManager.removeWaitCursor(token);
        }


        return "OK";
    }

    private void createSnapshot(String filename) {
        IGV mainFrame = IGV.getFirstInstance();

        if (filename == null) {
            String locus = FrameManager.getDefaultFrame().getFormattedLocusString();
            filename = locus.replaceAll(":", "_").replace("-", "_") + ".png";
        }

        File file = snapshotDirectory == null ? new File(filename) : new File(snapshotDirectory, filename);
        System.out.println("Snapshot: " + file.getAbsolutePath());

        SnapshotUtilities.doSnapshotOffscreen(mainFrame.getMainPanel(), file);
    }

    private static RegionScoreType getRegionSortOption(String str) {
        if (str == null) return null;
        String option = str.toUpperCase();
        try {
            return RegionScoreType.valueOf(option);
        }
        catch (Exception e) {
            return null;
        }
    }


    //START, STRAND, NUCLEOTIDE, QUALITY, SAMPLE, READ_GROUP
    private static AlignmentTrack.SortOption getAlignmentSortOption(String str) {
        String option = str.toLowerCase();
        if (str == null || option.equals("base")) {
            return AlignmentTrack.SortOption.NUCELOTIDE;
        } else if (option.equals("strand")) {
            return AlignmentTrack.SortOption.STRAND;

        } else if (option.equals("start") || option.equals("position")) {
            return AlignmentTrack.SortOption.START;

        } else if (option.equals("quality")) {
            return AlignmentTrack.SortOption.QUALITY;

        } else if (option.equals("sample")) {
            return AlignmentTrack.SortOption.SAMPLE;

        } else if (option.equals("readGroup") || option.equals("read_group")) {
            return AlignmentTrack.SortOption.READ_GROUP;
        } else if (option.equals("insertSize") || option.equals("insert_size")) {
            return AlignmentTrack.SortOption.INSERT_SIZE;
        }
        return AlignmentTrack.SortOption.NUCELOTIDE;
    }
}
