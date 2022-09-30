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

import org.broad.igv.Globals;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.oauth.OAuthUtils;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.SortOption;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionReader;
import org.broad.igv.session.SessionWriter;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.OverlayTracksMenuAction;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.SnapshotUtilities;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.*;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.net.URLDecoder;
import java.nio.charset.Charset;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;

public class CommandExecutor {

    private static Logger log = LogManager.getLogger(CommandExecutor.class);

    private File snapshotDirectory;
    private IGV igv;
    private String scriptDir;
    private int sleepInterval = 0; //2000;

    public CommandExecutor(IGV igv) {
        this(igv, null);
    }

    public CommandExecutor(IGV igv, String scriptDir) {
        this.igv = igv;
        this.scriptDir = scriptDir;
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


    public String execute(String commandLine) {

        List<String> args = getArgs(StringUtils.breakQuotedString(commandLine, ' ').toArray(new String[]{}));

        String result = "OK";

        log.debug("Executing: " + commandLine);
        try {
            if (args.size() == 0) {
                return "Empty command string";
            }

            String cmd = args.get(0).toLowerCase();
            String param1 = args.size() > 1 ? args.get(1) : null;
            String param2 = args.size() > 2 ? args.get(2) : null;
            String param3 = args.size() > 3 ? args.get(3) : null;
            String param4 = args.size() > 4 ? args.get(4) : null;

            if (cmd.equalsIgnoreCase("echo")) {
                result = param1 != null ? param1 : cmd;
            } else if (cmd.equalsIgnoreCase("gotoimmediate") || cmd.equalsIgnoreCase("goto")) {
                result = goto1(args);
            } else if (cmd.equalsIgnoreCase("addframes")) {
                result = addFrames(args);
            } else if (cmd.equalsIgnoreCase("scrolltotrack") || cmd.equalsIgnoreCase("gototrack")) {
                boolean res = this.igv.scrollToTrack(StringUtils.stripQuotes(param1));
                result = res ? "OK" : String.format("Error: Track %s not found", param1);
            } else if (cmd.equalsIgnoreCase("snapshotdirectory")) {
                result = setSnapshotDirectory(param1);
            } else if (cmd.equalsIgnoreCase("snapshot")) {
                result = createSnapshot(param1, param2);
            } else if (cmd.equalsIgnoreCase("savesession")) {
                String filename = param1;
                result = saveSession(filename);
            } else if ((cmd.equalsIgnoreCase("loadfile") || cmd.equalsIgnoreCase("load")) && param1 != null) {
                result = load(param1, args.size() > 2 ? args.subList(2, args.size()) : null);
            } else if (cmd.equalsIgnoreCase("genome") && args.size() > 1) {
                result = genome(param1);
            } else if (cmd.equalsIgnoreCase("new") || cmd.equalsIgnoreCase("reset") || cmd.equalsIgnoreCase("clear")) {
                igv.newSession();
            } else if (cmd.equalsIgnoreCase("region")) {
                defineRegion(param1, param2, param3, param4);
            } else if (cmd.equalsIgnoreCase("sort")) {
                result = sort(param1, param2, param3, param4);
            } else if (cmd.equalsIgnoreCase("group")) {
                result = group(param1, param2);
            } else if (cmd.equalsIgnoreCase("colorBy")) {
                result = colorBy(param1, param2);
            } else if (cmd.equalsIgnoreCase("collapse")) {
                String trackName = parseTrackName(param1);
                setTrackDisplayMode(Track.DisplayMode.COLLAPSED, trackName);
            } else if (cmd.equalsIgnoreCase("setSequenceStrand")) {
                setSequenceTrackStrand(Strand.fromString(param1));
            } else if (cmd.equalsIgnoreCase("setSequenceShowTranslation")) {
                boolean showTranslation;
                try {
                    if (param1.equalsIgnoreCase("true") || param1.equalsIgnoreCase("false")) {
                        showTranslation = Boolean.valueOf(param1);
                    } else {
                        return "ERROR: showTranslation value (" + param1 + ")is not 'true' or 'false'.";
                    }
                } catch (IllegalArgumentException e) {
                    return e.getMessage();
                }
                setSequenceShowTranslation(showTranslation);
            } else if (cmd.equalsIgnoreCase("expand")) {
                String trackName = parseTrackName(param1);
                setTrackDisplayMode(Track.DisplayMode.EXPANDED, trackName);
            } else if (cmd.equalsIgnoreCase("squish")) {
                String trackName = parseTrackName(param1);
                setTrackDisplayMode(Track.DisplayMode.SQUISHED, trackName);
            } else if (cmd.equalsIgnoreCase("remove")) {
                String trackName = parseTrackName(param1);
                result = removeTrack(trackName);
            } else if (cmd.equalsIgnoreCase("tweakdivider")) {
                igv.tweakPanelDivider();
            } else if (cmd.equalsIgnoreCase("setDataRange")) {
                result = this.setDataRange(param1, param2);
            } else if (cmd.equalsIgnoreCase("setLogScale")) {
                result = this.setLogScale(param1, param2);
            } else if (cmd.equalsIgnoreCase("setColor")) {
                result = this.setColor(param1, param2, false);
            } else if (cmd.equalsIgnoreCase("setAltColor")) {
                result = this.setColor(param1, param2, true);
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
            } else if ("oauth".equals(cmd)) {
                OAuthUtils.getInstance().getProvider().setAccessToken(param1);
            } else if (cmd.equalsIgnoreCase("sortByAttribute")) {
                result = sortByAttribute(args);
            } else if (cmd.equalsIgnoreCase("fitTracks")) {
                igv.fitTracksToPanel();
            } else if (cmd.equalsIgnoreCase("showAttributes")) {
                result = this.showAttributes(args);
            } else if (cmd.equalsIgnoreCase("showDataRange")) {
                result = this.setShowDataRange(param1, param2);
            } else if (cmd.equalsIgnoreCase("setTrackHeight")) {
                result = this.setTrackHeight(param1, param2);
            } else if (cmd.equalsIgnoreCase("overlay")) {
                result = this.overlay(args);
            } else if (cmd.equalsIgnoreCase("separate")) {
                result = this.separate(param1);
            } else {
                result = "UNKOWN COMMAND: " + commandLine;
                log.warn(result);
                return result;
            }

            igv.repaint();

            if (RuntimeUtils.getAvailableMemoryFraction() < 0.5) {
                log.debug("Running garbage collection");
                System.gc();
            }
            log.debug("Finished execution: " + commandLine + "  sleeping ....");
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

        log.debug(result);

        return result;
    }

    /**
     * Sort tracks by one or more sample attribute values.
     *
     * @param args
     * @return
     */
    private String sortByAttribute(List<String> args) {

        int nattributes = (args.size() - 1) / 2;
        if (nattributes == 0 || (args.size() - 1) % 2 != 0) {
            return String.format("Error: sortByAttribute usage: sortByAttribute attributeName asc|desc");
        }

        // Build a hash to support case insensitve attribute name comparison
        List<String> allAttributes = AttributeManager.getInstance().getAttributeNames();
        Map<String, String> attributeMap = new HashMap<>();
        for (String att : allAttributes) {
            attributeMap.put(att.toUpperCase(), att);
        }

        boolean[] ascending = new boolean[nattributes];
        String[] attributes = new String[nattributes];
        for (int attributeIndex = 0, i = 1; attributeIndex < nattributes; attributeIndex++, i += 2) {
            String attributeArg = StringUtils.stripQuotes(args.get(i)).toUpperCase();
            String attributeName = attributeMap.get(attributeArg);
            if (attributeName == null) {
                return String.format("Error: Attribute %s not found", attributeName);
            }
            String order = args.get(i + 1);
            ascending[attributeIndex] = order.equalsIgnoreCase("asc");
            attributes[attributeIndex] = attributeName;
        }
        igv.sortAllTracksByAttributes(attributes, ascending);
        return "OK";
    }

    private String removeTrack(String trackName) {
        if (trackName == null) {
            return "Error: NULL TRACK NAME";
        }
        List<Track> tracks = tracksMatchingName(trackName);
        if (tracks.size() > 0) {
            igv.deleteTracks(tracks);
            return "OK";
        }
        return String.format("Error: Track %s not found", trackName);
    }

    static String parseTrackName(String param1) {
        return param1 == null ? null : StringUtils.stripQuotes(param1);
    }

    private String overridePreference(String prefKey, String prefVal) {
        PreferencesManager.setOverride(prefKey, prefVal);
        return "OK";
    }

    private String setDataRange(String dataRangeString, String trackName) {

        String[] tokens = dataRangeString.split(",");
        //Min,max or min,baseline,max
        DataRange range = null;
        boolean autoscale =
                (dataRangeString.trim().equalsIgnoreCase("auto") ||
                        dataRangeString.trim().equalsIgnoreCase("autoscale"));

        if (!autoscale) {
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
        }
        List<Track> tracks = tracksMatchingName(trackName);
        for (Track track : tracks) {
            if (!autoscale) {
                if (track.getDataRange().isLog()) {
                    range.setType(DataRange.Type.LOG);  // Maintain "logness"
                }
                track.setDataRange(range);
            }
            track.setAutoScale(autoscale);

        }
        igv.repaint();
        return "OK";
    }

    private String setLogScale(String logScaleString, String trackName) {

        boolean logScale;
        try {
            if (logScaleString.equalsIgnoreCase("true") || logScaleString.equalsIgnoreCase("false")) {
                logScale = Boolean.valueOf(logScaleString);
            } else {
                return "ERROR: logscale value (" + logScaleString + ")is not 'true' or 'false'.";
            }
        } catch (IllegalArgumentException e) {
            return e.getMessage();
        }
        DataRange.Type scaleType = logScale == true ?
                DataRange.Type.LOG :
                DataRange.Type.LINEAR;
        List<Track> tracks = tracksMatchingName(trackName);
        for (Track track : tracks) {
            track.getDataRange().setType(scaleType);
        }
        igv.repaint(tracks);
        return "OK";
    }

    private String setColor(String colorString, String trackName, boolean alt) {

        try {
            Color color = ColorUtilities.stringToColorNoDefault(colorString);
            if (color == null) {
                return "Error: unrecognized color value " + colorString;
            }
            List<Track> tracks = tracksMatchingName(trackName);
            List<Track> affectedTracks = new ArrayList<>();
            for (Track track : tracks) {
                if (alt) {
                    track.setAltColor(color);
                } else {
                    track.setColor(color);
                }
                affectedTracks.add(track);

            }
            igv.repaint(affectedTracks);
            return "OK";
        } catch (Exception e) {
            return "Error setting track color: " + e.getMessage();
        }
    }

    private String overlay(List<String> args) {
        if (args.size() > 2) {
            String name = StringUtils.isQuoted(args.get(1)) ?
                    StringUtils.stripQuotes(args.get(1)) :
                    URLDecoder.decode(args.get(1), Charset.defaultCharset());
            List<DataTrack> tracks = new ArrayList<>();
            for (int i = 2; i < args.size(); i++) {
                List<Track> tmp = this.tracksMatchingName(args.get(i));
                for (Track t : tmp) {
                    if (t instanceof DataTrack) {
                        tracks.add((DataTrack) t);
                    }
                }
            }
            OverlayTracksMenuAction.merge(tracks, name);
            return "OK";
        } else {
            return "overlay command requires at least 2 arguments (trackName track1 track2 ...)";
        }
    }

    private String separate(String trackName) {
        List<Track> tmp = this.tracksMatchingName(trackName);
        if (tmp != null && tmp.size() > 0) {
            OverlayTracksMenuAction.unmerge(tmp);
            return "OK";
        } else {
            return "No track found matching " + trackName;
        }
    }

    private List<Track> tracksMatchingName(String name) {
        List<Track> tracks = igv.getAllTracks();
        if (name == null) {
            return tracks;
        } else {
            String altName = StringUtils.isQuoted(name) ?
                    StringUtils.stripQuotes(name) :
                    URLDecoder.decode(name, Charset.defaultCharset());

            return tracks.stream()
                    .filter(t ->
                            name.equalsIgnoreCase(t.getName()) || altName.equalsIgnoreCase(t.getName()) ||
                                    (t.getResourceLocator() != null && name.equals(t.getResourceLocator().getPath())) ||
                                    (t.getResourceLocator() != null && altName.equals(t.getResourceLocator().getPath()))
                    )
                    .collect(Collectors.toList());
        }
    }

    private String setViewAsPairs(String vAPString, String trackName) {
        List<Track> tracks = tracksMatchingName(trackName);
        boolean vAP = "false".equalsIgnoreCase(vAPString) ? false : true;
        for (Track track : tracks) {
            if (track instanceof AlignmentTrack) {
                AlignmentTrack atrack = (AlignmentTrack) track;
                atrack.setViewAsPairs(vAP);
            }
        }
        return "OK";
    }

    private String showAttributes(List<String> args) {
        // provides whitelist of visible attributes
        Set<String> hiddenAttributes = new HashSet<>(AttributeManager.getInstance().getAttributeNames());
        hiddenAttributes.addAll(igv.getSession().getHiddenAttributes());
        // Build a hash to support case insensitive attribute name comparison
        Map<String, String> attributeMap = new HashMap<>();
        for (String att : hiddenAttributes) {
            attributeMap.put(att.toUpperCase(), att);
        }
        for (int i = 1; i < args.size(); i++) {
            String attributeArg = StringUtils.stripQuotes(args.get(i)).toUpperCase();
            String attributeName = attributeMap.get(attributeArg);
            if (!hiddenAttributes.contains(attributeName)) {
                return String.format("Error: Attribute %s not found", attributeName);
            }
            hiddenAttributes.remove(attributeName);
        }
        igv.getSession().setHiddenAttributes(hiddenAttributes);
        igv.revalidateTrackPanels();
        return "OK";
    }


    private void setTrackDisplayMode(Track.DisplayMode mode, String trackName) {
        for (Track t : tracksMatchingName(trackName)) {
            t.setDisplayMode(mode);
        }
    }

    private void setSequenceTrackStrand(Strand trackStrand) {
        for (Track t : igv.getAllTracks()) {
            if (t instanceof SequenceTrack) {
                ((SequenceTrack) t).setStrand(trackStrand);
            }
        }
    }

    private void setSequenceShowTranslation(boolean shouldShowTranslation) {
        for (Track t : igv.getAllTracks()) {
            if (t instanceof SequenceTrack) {
                ((SequenceTrack) t).setShowTranslation(shouldShowTranslation);
            }
        }
    }

    private String setTrackHeight(String param1, String param2) {

        int height;
        String trackName;
        try {
            height = Integer.parseInt(param1);
            trackName = parseTrackName(param2);
        } catch (NumberFormatException e) {
            height = Integer.parseInt(param2);
            trackName = parseTrackName(param1);
        }

        height = Math.max(0, height);
        List<Track> tracks = tracksMatchingName(trackName);
        if (tracks.size() > 0) {
            for (Track track : tracks) {
                track.setHeight(height, true);
                igv.repaint(track);
            }
            return "OK";
        } else {
            return String.format("Error: Track %s not found", trackName);
        }
    }

    private String setShowDataRange(String show, String trackName) {

        boolean showDataRange;
        try {
            if (show.equalsIgnoreCase("true") || show.equalsIgnoreCase("false")) {
                showDataRange = Boolean.valueOf(show);
            } else {
                return "ERROR: showDataRange value (" + show + ") is not 'true' or 'false'.";
            }
        } catch (IllegalArgumentException e) {
            return e.getMessage();
        }
        List<Track> affectedTracks = new ArrayList<>();
        List<Track> tracks = tracksMatchingName(trackName);
        for (Track track : tracks) {
            if (track instanceof DataTrack) {
                ((DataTrack) track).setShowDataRange(showDataRange);
                affectedTracks.add(track);
            }
        }
        igv.repaint(affectedTracks);
        return "OK";
    }

    private String setSamplingWindowSize(String windowSize) {
        try {
            Integer.parseInt(windowSize);
            PreferencesManager.getPreferences().override(Constants.SAM_SAMPLING_WINDOW, String.valueOf(windowSize));
            return "OK";
        } catch (NumberFormatException e) {
            return "ERROR: SAMPLING WINDOW IS NOT A NUMBER: " + windowSize;
        }

    }

    private String setSamplingReadCount(String samplingReadCount) {
        try {
            Integer.parseInt(samplingReadCount);
            PreferencesManager.getPreferences().override(Constants.SAM_SAMPLING_COUNT, String.valueOf(samplingReadCount));
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

    public String setSleepInterval(String param1) {
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
        if (!ParsingUtils.fileExists(genomePath)) {
            String workingDirectory = System.getProperty("user.dir", "");
            genomePath = FileUtils.getAbsolutePath(genomeID, workingDirectory);
        }
        if (ParsingUtils.fileExists(genomePath)) {
            try {
                GenomeManager.getInstance().loadGenome(genomePath, null);
            } catch (IOException e) {
                result = "ERROR: Could not load genome: " + genomeID;
                MessageUtils.showMessage(result);
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
     * @param fileString path to file
     * @param params     optional parameters
     * @return
     * @throws IOException
     */
    private String load(String fileString, List<String> params) throws IOException {

        // remaining parameters might be "merge", "name", or "index"
        String name = null;
        String index = null;
        String coverage = null;
        String format = null;
        boolean merge = true;
        if (params != null) {
            for (String param : params) {
                if (param != null && param.startsWith("name=")) {
                    name = param.substring(5);
                } else if (param != null && param.startsWith("merge=")) {
                    String mergeString = param.substring(6);
                    merge = !mergeString.equalsIgnoreCase("false");
                } else if (param != null && param.startsWith("index=")) {
                    index = param.substring(6);
                } else if (param != null && param.startsWith("coverage=")) {
                    coverage = param.substring(9);
                } else if (param != null && param.startsWith("format=")) {
                    format = param.substring(7);
                }
            }
        }


        // Locus is not specified from port commands
        String locus = null;
        Map<String, String> ignore = null;
        String sort = null;
        String sortTag = null;
        return loadFiles(fileString, index, coverage, name, format, locus, merge, ignore, sort, sortTag, true);

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
     * @param dup
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
                     String sortTag,
                     boolean dup) {

        boolean isDataURL = ParsingUtils.isDataURL(fileString);

        List<String> files = isDataURL ? Arrays.asList(fileString) : breakFileString(fileString);
        List<String> indexFiles = isDataURL ? null : breakFileString(indexString);
        List<String> coverageFiles = breakFileString(coverageString);
        List<String> names = breakFileString(nameString);
        List<String> formats = breakFileString(formatString);
        if (names != null && names.size() != files.size()) {
            return "Error: If file is a comma-separated list, names must also be a comma-separated list of the same length";
        }
        if (indexFiles != null && indexFiles.size() != files.size()) {
            return "Error: If file is a comma-separated list, index must also be a comma-separated list of the same length";
        }
        if (isDataURL && formatString == null) {
            return "Error: format must be specified for dataURLs";
        }

        // Fix formats -- backward compatibility
        if (formats != null) {
            for (int i = 0; i < formats.size(); i++) {
                String formatOrExt = decodeFileString(formats.get(i));
                String format = formatOrExt.startsWith(".") ? formatOrExt.substring(1) : formatOrExt;
                formats.set(i, format);

            }
        }

        List<ResourceLocator> fileLocators = new ArrayList<>();

        if (!SessionReader.isSessionFile(fileString) && merge == false) {
            igv.newSession();
        }

        // Create set of loaded files
        Set<String> loadedFiles = new HashSet<>();
        for (ResourceLocator rl : igv.getDataResourceLocators()) {
            loadedFiles.add(rl.getPath());
        }


        // Loop through files
        for (int fi = 0; fi < files.size(); fi++) {
            String f = files.get(fi);
            if (!FileUtils.isRemote(f)) {
                File maybeFile = getFile(f);
                if (maybeFile.exists()) {
                    f = maybeFile.getAbsolutePath();
                } else {
                    maybeFile = new File(this.scriptDir, StringUtils.stripQuotes(f));
                    if (maybeFile.exists()) {
                        f = maybeFile.getAbsolutePath();
                    }
                }
            }

            if (isDataURL && formats == null) {
                return "Error: format must be specified for dataURLs";
            }


            // Skip already loaded files
            if (!dup && loadedFiles.contains(f)) continue;

            if (SessionReader.isSessionFile(f)) {
                igv.loadSession(f, locus);
            } else {

                ResourceLocator rl = new ResourceLocator(f);

                if (names != null) {
                    rl.setName(names.get(fi));
                } else if (isDataURL) {
                    rl.setName("Data");
                }
                if (indexFiles != null) {
                    rl.setIndexPath(indexFiles.get(fi));
                }
                if (coverageFiles != null) {
                    rl.setCoverage(coverageFiles.get(fi));
                }
                if (formats != null) {
                    rl.setFormat(formats.get(fi));
                }

                if (params != null) {
                    String trackLine = createTrackLine(params);
                    rl.setTrackLine(trackLine);
                }

                if (!isDataURL && rl.isLocal()) {
                    File file = new File(rl.getPath());
                    if (!file.exists()) {
                        return "Error: " + f + " does not exist.";
                    }
                }
                fileLocators.add(rl);
            }
        }


        if (fileLocators.size() > 0) {
            igv.loadTracks(fileLocators);
        }

        if (locus != null && !locus.equals("null")) {
            igv.goToLocus(locus);
            //If locus is a single base, we sort by base
            String[] tokens = locus.split(":", 2);
            if (tokens.length == 2) {
                try {
                    int pos = Integer.parseInt(tokens[1].replace(",", ""));
                    if (pos >= 0 && sort == null) sort = "base";
                } catch (Exception e) {
                    //pass
                }
            }
        }

        if (sort != null) {
            final SortOption sortOption = getAlignmentSortOption(sort);
            igv.sortAlignmentTracks(sortOption, sortTag, false);
        }

        return CommandListener.OK;
    }


    static List<String> breakFileString(String fileString) {
        if (fileString == null) {
            return null;
        }
        List<String> files = StringUtils.breakQuotedString(fileString, ',');
        if (files.size() == 1) {
            // String might be URL encoded
            List<String> files2 = StringUtils.breakQuotedString(fileString.replaceAll("%2C", ","), ',');
            if (files2.size() > 1) {
                files = files2;
            }
        }

        // URL decode strings.   This could be problematic as we don't know they were encoded, but its been like this
        // since 2009
        for (int ii = 0; ii < files.size(); ii++) {
            files.set(ii, decodeFileString(files.get(ii)));
        }

        return files;
    }


    /**
     * If {@code fileString} is a URL and can be decoded,
     * return the decoded version. Otherwise return the original.
     *
     * @param fileString
     * @return
     */
    static String decodeFileString(String fileString) {
        if (StringUtils.isQuoted(fileString)) {
            return StringUtils.stripQuotes(fileString);
        } else if (needsDecode(fileString)) {
            return StringUtils.decodeURL(fileString);
        } else {
            return fileString;
        }
    }

    static boolean needsDecode(String fileString) {

        String decodedString = decodeSafe(fileString);
        return (decodedString != null && (URLUtils.isURL(fileString) || URLUtils.isURL(decodedString)));
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

        File snapshotDir = getFile(param1);

        String result;
        if (snapshotDir.exists()) {
            snapshotDirectory = snapshotDir;
            result = "OK";
        } else {
            createParents(snapshotDir);
            snapshotDir.mkdir();
            if (snapshotDir.exists()) {
                snapshotDirectory = snapshotDir;
                result = "OK";
            } else {
                result = "ERROR: directory: " + param1 + " does not exist";
            }
        }
        return result;
    }

    /**
     * Fetch a file for the given path, which might be absolute, relative to the user home directory, or relative
     * to the script path
     *
     * @param path
     * @return
     */
    private File getFile(String path) {
        // Strip trailing & leading quotes
        path = StringUtils.stripQuotes(path);

        // Replace leading ~ with home directory
        if (path.startsWith("~")) {
            path = System.getProperty("user.home") + path.substring(1);
            return new File(path);
        } else if (path.startsWith("$SCRIPT_DIR") && this.scriptDir != null) {
            return new File(this.scriptDir + path.substring("$SCRIPT_DIR".length()));
        } else {
            File maybeFile = new File(path);
            if (maybeFile.exists()) {
                return maybeFile;
            } else if (scriptDir != null) {
                maybeFile = new File(this.scriptDir, path);
                if (maybeFile.exists()) {
                    return maybeFile;
                }
            }
        }
        return new File(path);
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

    private String addFrames(List<String> args) {
        if (args == null || args.size() < 3) {
            return "ERROR: missing locus parameter";
        }
        FrameManager.addFrames(args.subList(1, args.size()));
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


    private String sort(String sortArg, String param2, String param3, String param4) {

        RegionScoreType regionSortOption = getRegionSortOption(sortArg);
        if (regionSortOption != null) {
            // Segmented copy number
            RegionOfInterest roi = null;
            if (param2 != null) {
                Locus locus = Locus.fromString(param2);
                if (locus != null) {
                    int start = Math.max(0, locus.getStart() - 1);
                    roi = new RegionOfInterest(locus.getChr(), start, locus.getEnd(), "");
                }
            }
            igv.sortByRegionScore(roi, regionSortOption, FrameManager.getFirstFrame());
            return "OK";
        } else {
            // Alignments
            String tag = null;
            String locusString;
            String reverseString;
            if (sortArg.equalsIgnoreCase("tag")) {
                tag = param2;
                locusString = param3;
                reverseString = param4;
            } else {
                locusString = param2;
                reverseString =  param3;
            }

            // Special case, "reverse" is a resered word for inverting sorting.  Locus is optional
            if(reverseString == null && "reverse".equalsIgnoreCase(locusString)) {
                reverseString = locusString;
                locusString = null;
            }

            Double location = null;
            if (locusString != null && locusString.trim().length() > 0) {
                // Apparently there have been 2 conventions for "location", a full locus string and a base position
                // Try locus string first
                Locus locus = Locus.fromString(locusString);
                if (locus != null) {
                    location = (double) locus.getStart();
                } else {
                    try {
                        location = Double.valueOf(locusString.replace(",", ""));
                        if (location != null) location--;
                    } catch (NumberFormatException e) {
                        return "Error parsing location: " + locusString;
                    }
                }
            }

            boolean invertSort = "reverse".equalsIgnoreCase(reverseString);

            igv.sortAlignmentTracks(getAlignmentSortOption(sortArg), location, tag, invertSort);
            return "OK";
        }
    }


    private String group(String groupArg, String tagArg) {
        final AlignmentTrack.GroupOption groupOption = getAlignmentGroupOption(groupArg);
        Range r = null;
        if (groupOption == AlignmentTrack.GroupOption.BASE_AT_POS) {
            if (tagArg == null) {
                return "Error: position is required";
            } else {
                try {
                    Locus locus = Locus.fromString(tagArg);
                    r = new Range(locus.getChr(), locus.getStart(), locus.getEnd());
                } catch (Exception e) {
                    return "Error: invalid position: " + tagArg;
                }
            }
        }
        igv.groupAlignmentTracks(groupOption, tagArg, r);
        return "OK";
    }

    private String colorBy(String colorArg, String tagArg) {
        final AlignmentTrack.ColorOption colorOption = AlignmentTrack.ColorOption.valueOf(colorArg.toUpperCase());
        igv.colorAlignmentTracks(colorOption, tagArg);
        return "OK";
    }


    private String createSnapshot(String filename, String region) {

        if (filename == null) {
            String locus = FrameManager.getCurrentLocusString();
            filename = locus.replaceAll(":", "_").replace("-", "_") + ".png";
        } else {
            filename = StringUtils.stripQuotes(filename);
        }

        File file;
        if (snapshotDirectory == null) {
            file = getFile(filename);
            if (!file.getAbsoluteFile().getParentFile().exists()) {
                createParents(file);
            }
        } else {
            file = new File(snapshotDirectory, filename);
        }

        Component target = null;
        if (region == null || region.trim().length() == 0) {
            target = this.igv.getContentPane().getMainPanel();
        } else if ("trackpanels".equalsIgnoreCase(region)) {
            target = this.igv.getMainPanel().getCenterSplitPane();
        }

        if (target == null) {
            String msg = "ERROR. Could not create snapshot. Unknown region: " + region;
            log.error(msg);
            return msg;
        }

        try {
            return this.igv.createSnapshotNonInteractive(target, file, true);
        } catch (Exception e) {
            log.error(e);
            return e.getMessage();
        }
    }

    private String saveSession(String filename) {
        Session currentSession = igv.getSession();
        if (!filename.endsWith(".xml")) {
            filename = filename + ".xml";
        }
        File targetFile = getFile(filename);
        if (targetFile.getParentFile().exists()) {
            currentSession.setPath(targetFile.getAbsolutePath());
            try {
                (new SessionWriter()).saveSession(currentSession, targetFile);
                return "OK";
            } catch (Exception e) {
                return "Error writingin sesssion: " + e.getMessage();
            }
        } else {
            return "Error: directory not found: " + targetFile.getParentFile().getAbsolutePath();
        }
    }


    private static void createParents(File outputFile) {
        File parent = outputFile.getParentFile();
        if (parent != null && !parent.exists()) {
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

    //todo check that this matches all the existing sorts
    private static SortOption getAlignmentSortOption(String str) {
        str = str == null ? "base" : str;
        if (str.equalsIgnoreCase("start") || str.equalsIgnoreCase("position")) {
            return SortOption.START;
        } else if (str.equalsIgnoreCase("strand")) {
            return SortOption.STRAND;
        } else if (str.equalsIgnoreCase("base")) {
            return SortOption.NUCLEOTIDE;
        } else if (str.equalsIgnoreCase("quality")) {
            return SortOption.QUALITY;
        } else if (str.equalsIgnoreCase("sample")) {
            return SortOption.SAMPLE;
        } else if (str.equalsIgnoreCase("readGroup") || str.equalsIgnoreCase("read_group")) {
            return SortOption.READ_GROUP;
        } else if (str.equalsIgnoreCase("insertSize") || str.equalsIgnoreCase("insert_size")) {
            return SortOption.INSERT_SIZE;
        } else if (str.equalsIgnoreCase("firstOfPairStrand")) {
            return SortOption.FIRST_OF_PAIR_STRAND;
        } else if (str.equalsIgnoreCase("mateChr")) {
            return SortOption.MATE_CHR;
        } else if (str.equalsIgnoreCase("readOrder")) {
            return SortOption.READ_ORDER;
        } else if (str.equalsIgnoreCase("readname")) {
            return SortOption.READ_NAME;
        } else if (str.equalsIgnoreCase("alignedReadLength")) {
            return SortOption.ALIGNED_READ_LENGTH;
        } else {
            try {
                return SortOption.valueOf(str.toUpperCase());
            } catch (IllegalArgumentException e) {
                log.error("Unknown sort option: " + str);
                return SortOption.NUCLEOTIDE;
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
        } else if (str.equalsIgnoreCase("library")) {
            return AlignmentTrack.GroupOption.LIBRARY;
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


}
