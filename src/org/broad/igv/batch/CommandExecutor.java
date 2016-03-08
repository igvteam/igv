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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.batch;

import com.google.common.collect.Iterables;
import com.google.common.eventbus.EventBus;
import com.google.common.eventbus.Subscribe;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.dev.api.batch.Command;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.DataLoadedEvent;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.SnapshotUtilities;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.*;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

public class CommandExecutor {

    private static Logger log = Logger.getLogger(CommandExecutor.class);

    private File snapshotDirectory;
    private IGV igv;
    private int sleepInterval = 2000;


    public CommandExecutor() {
        igv = IGV.getInstance();
    }

    private List<String> getArgs(String[] tokens) {
        List<String> args = new ArrayList(tokens.length);
        for (String s : tokens) {
            if (s.trim().length() > 0) {
                args.add(s.trim());
            }
        }
        return args;
    }

    private String executeCustomCommand(String cmd, List<String> args) {

        List<String> subArgs = Collections.emptyList();
        if (args.size() > 1) subArgs = args.subList(1, args.size());
        try {
            Object ocmmand = RuntimeUtils.loadInstanceForName(cmd, null);
            Command command = (Command) ocmmand;
            return command.run(subArgs);
        } catch (ClassNotFoundException e) {
            return null;
        } catch (Exception e) {
            return e.getMessage();
        }
    }

    public String execute(String command) {

        List<String> args = getArgs(StringUtils.breakQuotedString(command, ' ').toArray(new String[]{}));

        String result = "OK";


        System.out.println();
        log.debug("Executing: " + command);
        try {
            if (args.size() == 0) {
                return "Empty command string";
            }
            //Custom command, user can make own
            String custRes = executeCustomCommand(args.get(0), args);
            if (custRes != null) {
                result = custRes;
            } else {

                String cmd = args.get(0).toLowerCase();
                String param1 = args.size() > 1 ? args.get(1) : null;
                String param2 = args.size() > 2 ? args.get(2) : null;
                String param3 = args.size() > 3 ? args.get(3) : null;
                String param4 = args.size() > 4 ? args.get(4) : null;


                if (cmd.equalsIgnoreCase("echo")) {
                    result = cmd;
                } else if (cmd.equalsIgnoreCase("gotoimmediate")) {
                    return gotoImmediate(args);
                } else if (cmd.equalsIgnoreCase("goto")) {
                    result = goto1(args);
                } else if (cmd.equalsIgnoreCase("gototrack")) {
                    boolean res = IGV.getInstance().scrollToTrack(param1);
                    result = res ? "OK" : String.format("Error: Track %s not found", param1);
                } else if (cmd.equalsIgnoreCase("snapshotdirectory")) {
                    result = setSnapshotDirectory(param1);
                } else if (cmd.equalsIgnoreCase("snapshot")) {
                    String filename = param1;
                    result = createSnapshot(filename, param2);
                } else if ((cmd.equalsIgnoreCase("loadfile") || cmd.equalsIgnoreCase("load")) && param1 != null) {
                    result = load(param1, param2, param3, param4);
                } else if (cmd.equalsIgnoreCase("genome") && args.size() > 1) {
                    result = genome(param1);
                } else if (cmd.equalsIgnoreCase("new") || cmd.equalsIgnoreCase("reset") || cmd.equalsIgnoreCase("clear")) {
                    igv.newSession();
                } else if (cmd.equalsIgnoreCase("region")) {
                    defineRegion(param1, param2, param3, param4);
                } else if (cmd.equalsIgnoreCase("sort")) {
                    sort(param1, param2, param3, param4);
                } else if (cmd.equalsIgnoreCase("group")) {
                    group(param1, param2);
                } else if (cmd.equalsIgnoreCase("collapse")) {
                    String trackName = parseTrackName(param1);
                    igv.setTrackDisplayMode(Track.DisplayMode.COLLAPSED, trackName);
                } else if (cmd.equalsIgnoreCase("expand")) {
                    String trackName = parseTrackName(param1);
                    igv.setTrackDisplayMode(Track.DisplayMode.EXPANDED, trackName);
                } else if (cmd.equalsIgnoreCase("squish")) {
                    String trackName = parseTrackName(param1);
                    igv.setTrackDisplayMode(Track.DisplayMode.SQUISHED, trackName);
                } else if (cmd.equalsIgnoreCase("remove")) {
                    String trackName = parseTrackName(param1);
                    result = removeTrack(trackName);
                } else if (cmd.equalsIgnoreCase("tweakdivider")) {
                    igv.tweakPanelDivider();
                } else if (cmd.equalsIgnoreCase("setDataRange")) {
                    result = this.setDataRange(param1, param2);
                } else if (cmd.equalsIgnoreCase("maxpanelheight") && param1 != null) {
                    return setMaxPanelHeight(param1);
                } else if (cmd.equalsIgnoreCase("tofront")) {
                    return UIUtilities.bringToFront();
                } else if (cmd.equalsIgnoreCase("viewaspairs")) {
                    return setViewAsPairs(param1, param2);
                } else if (cmd.equalsIgnoreCase("samplingwindowsize")) {
                    return this.setSamplingWindowSize(param1);
                } else if (cmd.equalsIgnoreCase("maxdepth") || (cmd.equalsIgnoreCase("samplingreadcount"))) {
                    return this.setSamplingReadCount(param1);
                } else if (cmd.equalsIgnoreCase("setSleepInterval")) {
                    return this.setSleepInterval(param1);
                } else if (cmd.equalsIgnoreCase("setCredentials")) {
                    return this.setCredentials(param1, param2);
                } else if (cmd.equalsIgnoreCase("clearCredentials")) {
                    return this.clearCredentials();
                } else if (cmd.equals("preference")) {
                    return this.overridePreference(param1, param2);
                } else if (cmd.equalsIgnoreCase("version")) {
                    return Globals.VERSION;
                } else if (cmd.equals("exit")) {
                    System.exit(0);
                } else if (cmd.equals("zoomin")) {
                    FrameManager.incrementZoom(1);
                } else if (cmd.equals("zoomout")) {
                    FrameManager.incrementZoom(-1);
                } else {
                    result = "UNKOWN COMMAND: " + command;
                    log.error(result);
                    return result;
                }
            }
            igv.doRefresh();

            if (RuntimeUtils.getAvailableMemoryFraction() < 0.5) {
                log.debug("Running garbage collection");
                System.gc();
            }
            log.debug("Finished execution: " + command + "  sleeping ....");
            if (sleepInterval > 0) try {
                Thread.sleep(sleepInterval);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            log.debug("Finished sleeping");

        } catch (IOException e) {
            log.error(e);
            result = "Error: " + e.getMessage();
        }
        log.info(result);

        return result;
    }

    private String removeTrack(String trackName) {
        if (trackName == null) return "Error: NULL TRACK NAME";
        for (Track track : igv.getAllTracks()) {
            if (track.getName().equals(trackName)) {
                igv.removeTracks(Arrays.asList(track));
                return "OK";
            }
        }
        return String.format("Error: Track %s not found", trackName);
    }

    static String parseTrackName(String param1) {
        return param1 == null ? null : StringUtils.stripQuotes(param1);
    }

    private String overridePreference(String prefKey, String prefVal) {
        PreferenceManager.getInstance().override(prefKey, prefVal);
        return "OK";
    }

    private String setDataRange(String dataRangeString, String trackName) {
        List<Track> tracks = igv.getAllTracks();
        String[] tokens = dataRangeString.split(",");
        //Min,max or min,baseline,max
        DataRange range;
        try {
            if (tokens.length == 2) {
                range = new DataRange(Float.parseFloat(tokens[0]), Float.parseFloat(tokens[1]));
            } else if (tokens.length == 3) {
                range = new DataRange(Float.parseFloat(tokens[0]), Float.parseFloat(tokens[1]), Float.parseFloat(tokens[2]));
            } else {
                throw new IllegalArgumentException(String.format("ERROR: parsing %s for data range. \n" +
                        "String must be of form <min,max> or <min,baseline,max>", dataRangeString));
            }
        } catch (NumberFormatException e) {
            return "ERROR: Could not parse input string as a Float. " + e.getMessage();
        } catch (IllegalArgumentException e) {
            return e.getMessage();
        }
        for (Track track : tracks) {
            if (trackName == null || trackName.equalsIgnoreCase(track.getName())) {
                track.setDataRange(range);
                track.setAutoScale(false);
            }
        }
        return "OK";
    }

    private String setViewAsPairs(String vAPString, String trackName) {
        List<Track> tracks = igv.getAllTracks();
        boolean vAP = "false".equalsIgnoreCase(vAPString) ? false : true;
        for (Track track : tracks) {
            if (track instanceof AlignmentTrack) {
                if (trackName == null || trackName.equalsIgnoreCase(track.getName())) {
                    AlignmentTrack atrack = (AlignmentTrack) track;
                    atrack.setViewAsPairs(vAP);
                }
            }
        }
        return "OK";
    }

    private String setSamplingWindowSize(String windowSize) {
        try {
            Integer.parseInt(windowSize);
            PreferenceManager.getInstance().override(PreferenceManager.SAM_SAMPLING_WINDOW, String.valueOf(windowSize));
            return "OK";
        } catch (NumberFormatException e) {
            return "ERROR: SAMPLING WINDOW IS NOT A NUMBER: " + windowSize;
        }

    }

    private String setSamplingReadCount(String samplingReadCount) {
        try {
            Integer.parseInt(samplingReadCount);
            PreferenceManager.getInstance().override(PreferenceManager.SAM_SAMPLING_COUNT, String.valueOf(samplingReadCount));
            return "OK";
        } catch (NumberFormatException e) {
            return "ERROR: SAMPLING READ COUNT IS NOT A NUMBER: " + samplingReadCount;
        }

    }

    private String gotoImmediate(List<String> args) {
        return goto1(args);
    }

    private String setMaxPanelHeight(String param1) {
        try {
            Integer h = Integer.parseInt(param1.trim());
            SnapshotUtilities.setMaxPanelHeight(h);
            return "OK";
        } catch (NumberFormatException e) {
            return "ERROR - max panel height value ('" + param1 + ".) must be an integer number";
        }
    }

    private String setSleepInterval(String param1) {
        try {
            sleepInterval = Integer.parseInt(param1.trim());
            return "OK";
        } catch (NumberFormatException e) {
            return "ERROR - sleep interval value ('" + param1 + ".) must be an integer number";
        }
    }

    private String setCredentials(String param1, String param2) {
        HttpUtils.getInstance().setDefaultUserName(param1);
        HttpUtils.getInstance().setDefaultPassword(param2);
        return "OK";
    }

    private String clearCredentials() {
        HttpUtils.getInstance().clearDefaultCredentials();
        return "OK";
    }


    private String genome(String param1) {
        if (param1 == null) {
            return "ERROR missing genome parameter";
        }
        String result = "OK";
        String genomeID = param1;

        igv.selectGenomeFromList(genomeID);
        if (GenomeManager.getInstance().getCurrentGenome().getId().equals(genomeID)) {
            return result;
        }


        String genomePath = genomeID;
        if (!ParsingUtils.pathExists(genomePath)) {
            String workingDirectory = System.getProperty("user.dir", "");
            genomePath = FileUtils.getAbsolutePath(genomeID, workingDirectory);
        }
        if (ParsingUtils.pathExists(genomePath)) {
            try {
                igv.loadGenome(genomePath, null, true);
            } catch (IOException e) {
                throw new RuntimeException("Error loading genome: " + genomeID);
            }
        } else {
            result = "ERROR: Could not locate genome: " + genomeID;
            MessageUtils.showMessage(result);
        }


        return result;
    }

    /**
     * Load function for port and batch script
     *
     * @param fileString
     * @param param2
     * @param param3
     * @return
     * @throws IOException
     */
    private String load(String fileString, String param2, String param3, String param4) throws IOException {

        // Default for merge is "true" for session files,  "false" otherwise
        String tmpFile = StringUtils.stripQuotes(fileString);
        boolean merge = !(tmpFile.endsWith(".xml") || tmpFile.endsWith(".php") || tmpFile.endsWith(".php3"));

        // remaining parameters might be "merge", "name", or "index"
        String name = null;
        String index = null;
        String coverage = null;
        String format = null;
        for (String param : Arrays.asList(param2, param3)) {
            if (param != null && param.startsWith("name=")) {
                name = param.substring(5);
            } else if (param != null && param.startsWith("merge=")) {
                String mergeString = param.substring(6);
                merge = mergeString.equalsIgnoreCase("true");
            } else if (param != null && param.startsWith("index=")) {
                index = param.substring(6);
            } else if (param != null && param.startsWith("coverage=")) {
                coverage = param.substring(9);
            } else if (param != null && param.startsWith("format=")) {
                format = param.substring(7);
            }
        }
        // Locus is not specified from port commands
        String locus = null;
        Map<String, String> params = null;
        return loadFiles(fileString, index, coverage, name, format, locus, merge, params);
    }

    String loadFiles(final String fileString,
                     final String indexString,
                     final String coverageString,
                     final String nameString,
                     final String formatString,
                     final String locus,
                     final boolean merge,
                     Map<String, String> params) throws IOException {
        return loadFiles(fileString, indexString, coverageString, nameString, formatString, locus, merge, params, null, null);
    }

    /**
     * Load files -- used by port, batch, and http commands
     *
     * @param fileString
     * @param locus
     * @param merge
     * @param nameString
     * @param params
     * @param sort
     * @param sortTag    Used iff sort == SortOption.TAG
     * @return
     * @throws IOException
     */
    String loadFiles(final String fileString,
                     final String indexString,
                     final String coverageString,
                     final String nameString,
                     final String formatString,
                     final String locus,
                     final boolean merge,
                     Map<String, String> params,
                     String sort,
                     String sortTag) throws IOException {


        log.debug("Run load files");

        List<String> files = StringUtils.breakQuotedString(fileString, ',');
        List<String> names = StringUtils.breakQuotedString(nameString, ',');
        List<String> indexFiles = StringUtils.breakQuotedString(indexString, ',');
        List<String> coverageFiles = StringUtils.breakQuotedString(coverageString, ',');
        List<String> formats = StringUtils.breakQuotedString(formatString, ',');

        if (files.size() == 1) {
            // String might be URL encoded
            files = StringUtils.breakQuotedString(fileString.replaceAll("%2C", ","), ',');
            names = nameString != null ? StringUtils.breakQuotedString(nameString.replaceAll("%2C", ","), ',') : null;
            indexFiles = indexString != null ? StringUtils.breakQuotedString(indexString.replaceAll("%2C", ","), ',') : null;
            coverageFiles = coverageString != null ? StringUtils.breakQuotedString(coverageString.replaceAll("%2C", ","), ',') : null;
            formats = formatString != null ? StringUtils.breakQuotedString(formatString.replaceAll("%2C", ","), ',') : null;
        }

        if (names != null && names.size() != files.size()) {
            return "Error: If file is a comma-separated list, names must also be a comma-separated list of the same length";
        }

        if (indexFiles != null && indexFiles.size() != files.size()) {
            return "Error: If file is a comma-separated list, index must also be a comma-separated list of the same length";
        }


        // Must decode URLs (local or remote), but leave local file paths only
        for (int ii = 0; ii < files.size(); ii++) {
            files.set(ii, decodeFileString(files.get(ii).replace("\"", "")));
            if (names != null) {
                names.set(ii, names.get(ii).replace("\"", ""));
            }
            if (indexFiles != null) {
                indexFiles.set(ii, decodeFileString(indexFiles.get(ii).replace("\"", "")));
            }
            if (coverageFiles != null) {
                coverageFiles.set(ii, decodeFileString(coverageFiles.get(ii).replace("\"", "")));
            }
            if (formatString != null) {
                formats.set(ii, decodeFileString(formats.get(ii).replace("\"", "")));
            }
        }

        List<ResourceLocator> fileLocators = new ArrayList<ResourceLocator>();
        List<String> sessionPaths = new ArrayList<String>();

        if (!merge) {
            // If this is a session file start fresh without asking, otherwise ask
            boolean unload = !merge;
            if (fileString.endsWith(".xml") || fileString.endsWith(".php") || fileString.endsWith(".php3")) {
                unload = !merge;
            } else {
                unload = MessageUtils.confirm("Unload current session before loading new tracks?");
            }
            if (unload) {
                igv.newSession();
            }
        }

        // Create set of loaded files
        Set<String> loadedFiles = new HashSet<String>();
        for (ResourceLocator rl : igv.getDataResourceLocators()) {
            loadedFiles.add(rl.getPath());
        }

        // Loop through files

        for (int fi = 0; fi < files.size(); fi++) {

            String f = files.get(fi);

            // Skip already loaded files TODO -- make this optional?  Check for change?
            if (loadedFiles.contains(f)) continue;

            if (f.endsWith(".xml") || f.endsWith(".php") || f.endsWith(".php3") || f.endsWith(".session")) {
                sessionPaths.add(f);
            } else {
                ResourceLocator rl = new ResourceLocator(f);

                if (rl.isLocal()) {
                    File file = new File(rl.getPath());
                    if (!file.exists()) {
                        return "Error: " + f + " does not exist.";
                    }
                }

                if (names != null) {
                    rl.setName(names.get(fi));
                }
                if (indexFiles != null) {
                    rl.setIndexPath(indexFiles.get(fi));
                }
                if (coverageFiles != null) {
                    rl.setCoverage(coverageFiles.get(fi));
                }
                if (formats != null) {
                    String format = formats.get(fi);
                    if (!format.startsWith(".")) format = "." + format;
                    rl.setType(format);
                }
                if (params != null) {
                    String trackLine = createTrackLine(params);
                    rl.setTrackLine(trackLine);
                }
                fileLocators.add(rl);
            }
        }

        for (String sessionPath : sessionPaths) {
            igv.restoreSessionSynchronous(sessionPath, locus, merge);
        }

        final Future loadTask = igv.loadTracks(fileLocators);

        if (locus != null && !locus.equals("null")) {
            igv.goToLocus(locus);
            //If locus is a single base, we sort by base
            String[] tokens = locus.split(":", 2);
            if (tokens.length == 2) {
                String chr = tokens[0];
                try {
                    int pos = Integer.parseInt(tokens[1].replace(",", ""));
                    if (pos >= 0 && sort == null) sort = "base";

                } catch (Exception e) {
                    //pass
                }
            }

        }

        if (sort != null) {
            submitPerformSort(loadTask, sort, sortTag);
        }

        return CommandListener.OK;
    }

    private void submitPerformSort(final Future loadTask, final String sort, final String sortTag) {
        final AlignmentTrack.SortOption sortOption = getAlignmentSortOption(sort);
        Runnable runnable = new Runnable() {
            @Override
            public void run() {
                try {
                    //We need to wait until the track is loaded. If loadTask is null,
                    //it was loaded synchronously
                    if (loadTask != null) {
                        Object res = loadTask.get();
                    }
                    //Thought we were done waiting, huh? Guess again
                    //Alignment tracks load alignment data asynchronously from the track

                    //Since sorting applies to all tracks, we only need to have 1 handler
                    AlignmentTrack track = null;
                    try {
                        track = Iterables.filter(igv.getAllTracks(), AlignmentTrack.class).iterator().next();
                        EventBus bus = track.getDataManager().getEventBus();
                        bus.register(new SortAlignmentsHandler(igv, bus, sortOption, sortTag));
                    } catch (NoSuchElementException e) {
                        //No alignment tracks found.
                        log.warn("Sort argument provided but no alignment tracks found");
                    }
                } catch (InterruptedException e) {
                    log.error(e.getMessage(), e);
                } catch (ExecutionException e) {
                    log.error(e.getMessage(), e);
                }

            }
        };
        LongRunningTask.submit(runnable);
    }

    /**
     * If {@code fileString} is a URL and can be decoded,
     * return the decoded version. Otherwise return the original.
     *
     * @param fileString
     * @return
     */
    static String decodeFileString(String fileString) {
        if (needsDecode(fileString)) {
            return StringUtils.decodeURL(fileString);
        } else {
            return fileString;
        }
    }

    static boolean needsDecode(String fileString) {
        String decodedString = decodeSafe(fileString);
        return (decodedString != null && (HttpUtils.isURL(fileString) || HttpUtils.isURL(decodedString)));
    }

    private static String decodeSafe(String string) {
        String tmp = null;
        try {
            tmp = StringUtils.decodeURL(string);
        } catch (Exception e) {
            log.warn(e.getMessage());
        }
        return tmp;
    }

    /**
     * Convert the parameter map to a UCSC track line.
     *
     * @param params
     * @return
     */
    private String createTrackLine(Map<String, String> params) {
        return params.get("hgt.customText");
//        StringBuffer buf = new StringBuffer();
//        buf.append("track ");
//        for(Map.Entry<String, String> entry : params.entrySet()) {
//            buf.append(entry.getKey());
//            buf.append("=");
//            buf.append(entry.getValue());
//            buf.append(" ");
//        }
//        return buf.toString();
    }


    /**
     * Set a directory to deposit image snapshots
     *
     * @param param1
     * @return
     */
    String setSnapshotDirectory(String param1) {

        if (param1 == null) {
            return "ERROR: missing directory parameter";
        }

        param1 = StringUtils.stripQuotes(param1);

        File parentDir = null;
        try {
            parentDir = getFile(param1);
        } catch (URISyntaxException e) {
            log.error("Error parsing directory path: " + param1, e);
            return "Error parsing directory path: " + param1;
        }

        String result;
        if (parentDir.exists()) {
            snapshotDirectory = parentDir;
            result = "OK";
        } else {
            createParents(parentDir);
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

    private File getFile(String param1) throws URISyntaxException {

        // Strip trailing & leading quotes
        if (param1.startsWith("\"")) param1 = param1.substring(1);
        if (param1.endsWith("\"")) param1 = param1.substring(0, param1.lastIndexOf('"'));

        // See if file contains spaces, if not no special treatment is required
        if (param1.indexOf(' ') < 0) {
            return new File(param1);
        } else {
            // If file is absolute use a URI,
            File f = new File(param1);
            if (f.isAbsolute()) {
                URI outputURI = new URI(("file://" + param1.replaceAll(" ", "%20")));
                return new File(outputURI);
            } else {
                return f;
            }
        }
    }

    private String goto1(List<String> args) {
        if (args == null || args.size() < 2) {
            return "ERROR: missing locus parameter";
        }
        String locus = args.get(1);
        for (int i = 2; i < args.size(); i++) {
            locus += (" " + args.get(i));
        }
        igv.goToLocus(locus);
        return "OK";
    }


    private void defineRegion(String param1, String param2, String param3, String param4) {

        RegionOfInterest roi = null;
        if (param1 != null && param2 != null && param3 != null) {
            int start = Math.max(0, Integer.parseInt(param2) - 1);
            int end = Integer.parseInt(param3);
            String desc = param4 != null ? param4 : "";
            roi = new RegionOfInterest(param1, start, end, desc);
        }
        if (param1 != null) {
            Locus locus = Locus.fromString(param1);
            if (locus != null) {
                int start = Math.max(0, locus.getStart() - 1);
                String desc = param2 != null ? param2 : "";
                roi = new RegionOfInterest(locus.getChr(), start, locus.getEnd(), desc);

            }
        }
        if (roi != null) {
            igv.addRegionOfInterest(roi);
        }
    }


    private void sort(String sortArg, String locusString, String param3, String param4) {

        RegionScoreType regionSortOption = getRegionSortOption(sortArg);
        String tag = "";
        if (regionSortOption != null) {
            RegionOfInterest roi = null;
            if (locusString != null) {
                Locus locus = Locus.fromString(locusString);
                if (locus != null) {
                    int start = Math.max(0, locus.getStart() - 1);
                    roi = new RegionOfInterest(locus.getChr(), start, locus.getEnd(), "");
                }
            }
            igv.sortByRegionScore(roi, regionSortOption, FrameManager.getDefaultFrame());

        } else {
            Double location = null;
            if (param3 != null && param3.trim().length() > 0) {
                try {
                    location = new Double(param3.replace(",", ""));
                    tag = param4;
                } catch (NumberFormatException e) {
                    tag = param3;
                }
            } else if (locusString != null && locusString.trim().length() > 0) {
                try {
                    location = new Double(locusString.replace(",", ""));
                    tag = param4;
                } catch (NumberFormatException e) {
                    tag = param3;
                }
            }

            //Convert from 1-based to 0-based
            if (location != null) location--;

            igv.sortAlignmentTracks(getAlignmentSortOption(sortArg), location, tag);

        }
        igv.repaintDataPanels();
    }

    private void group(String groupArg, String tagArg) {
        igv.groupAlignmentTracks(getAlignmentGroupOption(groupArg), tagArg);
        igv.repaintDataPanels();
    }


    private String createSnapshot(String filename, String region) {

        if (filename == null) {
            String locus = FrameManager.getDefaultFrame().getFormattedLocusString();
            filename = locus.replaceAll(":", "_").replace("-", "_") + ".png";
        } else {
            filename = StringUtils.stripQuotes(filename);
        }

        File file;
        if (snapshotDirectory == null) {
            try {
                file = getFile(filename);
                if (!file.getAbsoluteFile().getParentFile().exists()) {
                    createParents(file);
                }
            } catch (URISyntaxException e) {
                log.error("Error parsing directory path: " + filename, e);
                return "Error parsing directory path: " + filename;
            }
        } else {
            file = new File(snapshotDirectory, filename);
        }
        System.out.println("Snapshot: " + file.getAbsolutePath());


        Component target = null;
        if (region == null || region.trim().length() == 0) {
            target = IGV.getInstance().getContentPane().getMainPanel();
        } else if ("trackpanels".equalsIgnoreCase(region)) {
            target = IGV.getInstance().getMainPanel().getCenterSplitPane();
        }

        if (target == null) {
            String msg = "ERROR. Could not create snapshot. Unknown region: " + region;
            log.error(msg);
            return msg;
        }

        try {
            return IGV.getInstance().createSnapshotNonInteractive(target, file, true);
        } catch (IOException e) {
            log.error(e);
            return e.getMessage();
        }
    }


    private static void createParents(File outputFile) {
        File parent = outputFile.getParentFile();
        if (!parent.exists()) {
            parent.mkdirs();
        }

    }

    private static RegionScoreType getRegionSortOption(String str) {
        if (str == null) return null;
        String option = str.toUpperCase();
        try {
            return RegionScoreType.valueOf(option);
        } catch (IllegalArgumentException e) {
            return null;
        }
    }

    private static AlignmentTrack.SortOption getAlignmentSortOption(String str) {
        str = str == null ? "base" : str;
        if (str.equalsIgnoreCase("start") || str.equalsIgnoreCase("position")) {
            return AlignmentTrack.SortOption.START;
        } else if (str.equalsIgnoreCase("strand")) {
            return AlignmentTrack.SortOption.STRAND;
        } else if (str.equalsIgnoreCase("base")) {
            return AlignmentTrack.SortOption.NUCLEOTIDE;
        } else if (str.equalsIgnoreCase("quality")) {
            return AlignmentTrack.SortOption.QUALITY;
        } else if (str.equalsIgnoreCase("sample")) {
            return AlignmentTrack.SortOption.SAMPLE;
        } else if (str.equalsIgnoreCase("readGroup") || str.equalsIgnoreCase("read_group")) {
            return AlignmentTrack.SortOption.READ_GROUP;
        } else if (str.equalsIgnoreCase("insertSize") || str.equalsIgnoreCase("insert_size")) {
            return AlignmentTrack.SortOption.INSERT_SIZE;
        } else if (str.equalsIgnoreCase("firstOfPairStrand")) {
            return AlignmentTrack.SortOption.FIRST_OF_PAIR_STRAND;
        } else if (str.equalsIgnoreCase("mateChr")) {
            return AlignmentTrack.SortOption.MATE_CHR;
        } else {
            try {
                return AlignmentTrack.SortOption.valueOf(str.toUpperCase());
            } catch (IllegalArgumentException e) {
                log.error("Unknown sort option: " + str);
                return AlignmentTrack.SortOption.NUCLEOTIDE;
            }
        }

    }

    //      STRAND, SAMPLE, READ_GROUP, FIRST_OF_PAIR_STRAND, TAG, PAIR_ORIENTATION, MATE_CHROMOSOME, NONE, SUPPLEMENTARY

    private static AlignmentTrack.GroupOption getAlignmentGroupOption(String str) {
        if (str == null || str.length() == 0) {
            return AlignmentTrack.GroupOption.NONE;
        } else if (str.equalsIgnoreCase("strand")) {
            return AlignmentTrack.GroupOption.STRAND;
        } else if (str.equalsIgnoreCase("sample")) {
            return AlignmentTrack.GroupOption.SAMPLE;
        } else if (str.equalsIgnoreCase("readGroup") || str.equalsIgnoreCase("read_group")) {
            return AlignmentTrack.GroupOption.READ_GROUP;
        } else {
            try {
                return AlignmentTrack.GroupOption.valueOf(str.toUpperCase());
            } catch (IllegalArgumentException e) {
                log.error("Unknown group by option: " + str);
                return AlignmentTrack.GroupOption.NONE;
            }

        }

    }

    private static class SortAlignmentsHandler {

        private IGV igv = null;
        private EventBus bus = null;
        private AlignmentTrack.SortOption sortOption;
        private String sortTag;

        SortAlignmentsHandler(IGV igv, EventBus bus, AlignmentTrack.SortOption sortOption, String sortTag) {
            this.igv = igv;
            this.bus = bus;
            this.sortOption = sortOption;
            this.sortTag = sortTag;
        }

        @Subscribe
        public void received(DataLoadedEvent event) {
            boolean sorted = igv.sortAlignmentTracks(sortOption, sortTag);
            if (sorted) this.bus.unregister(this);
        }

    }
}
