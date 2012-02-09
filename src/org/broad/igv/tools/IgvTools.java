/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
package org.broad.igv.tools;


import jargs.gnu.CmdLineParser;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.RollingFileAppender;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.AbstractFeatureParser;
import org.broad.igv.feature.FeatureParser;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDescriptor;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.sam.reader.AlignmentIndexer;
import org.broad.igv.tdf.TDFUtils;
import org.broad.igv.tools.converters.ExpressionFormatter;
import org.broad.igv.tools.converters.GCTtoIGVConverter;
import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.converters.DensitiesToBedGraph;
import org.broad.igv.variant.util.VCFtoBed;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;

/**
 * Command line parser for "igvtools".
 *
 * @author jrobinso
 */
public class IgvTools {

    static String version = "@VERSION";

    static final String CMD_TILE = "tile";
    static final String CMD_TOTDF = "totdf";
    static final String CMD_COUNT = "count";
    static final String CMD_SORT = "sort";
    static final String CMD_INDEX = "index";
    static final String CMD_FORMATEXP = "formatexp";
    static final String CMD_VERSION = "version";
    static final String CMD_GUI = "gui";
    static final String CMD_HELP = "help";

    static Map<String, String> commandList = new HashMap<String, String>(9);

    static {
        String typ_end = "[options] inputFile outputFile genome\n";

        String tstring = "Usage: igvtools totdf " + typ_end;
        String zoom = "-z, --maxZoom \tMaximum zoom levels to precompute. Default 7\n";
        String wfunc = "-f, --windowFunctions \tComma-separated list of window functions. Possible values are min, mean, max\n";
        String probFile = "-p, --probeFile \tSpecifies a \"bed\" file to be used to map probe identifiers to locations.\n";
        tstring += zoom += wfunc += probFile;
        commandList.put(CMD_TOTDF, tstring);
        tstring = "Note: This command is deprecated and being renamed to totdf\n" + tstring;
        commandList.put(CMD_TILE, tstring);

        String countstring = "Usage: igvtools count " + typ_end;
        String wsize = "-w, --windowSize \tThe window size over which coverage is averaged. Defaults to 25 bp.\n";
        String ext = "-e, --extFactor \tThe read or feature is extended by the specified distance in bp prior to counting\n";
        //TODO MORE OPTIONS
        countstring += zoom += wfunc += wsize += ext;
        commandList.put(CMD_COUNT, countstring);

        String indexstring = "Usage: igvtools index inputFile\n";
        indexstring += "Index an alignment file.\n";
        indexstring += "Supported input file formats are: .sam, .aligned, .vcf, .psl, and .bed.\n";
        commandList.put(CMD_INDEX, indexstring);


        String sortstring = "Usage: igvtools sort [options] inputFile outputFile";
        String tmpDir = "-t, --tmpdir \tSpecify a temporary working directory. Default is users temp directory\n";
        String maxRec = "-m, --maxRecords \tThe maximum number of records to keep in memory during the sort\n";
        sortstring += tmpDir += maxRec;
        commandList.put(CMD_SORT, sortstring);

        String fexpstring = "Usage: igvtools formatexp inputFile outputFile\n";
        fexpstring += "Format GCT or RES files for display";
        commandList.put(CMD_FORMATEXP, fexpstring);

        commandList.put(CMD_VERSION, "Print igvtools version number");

        commandList.put(CMD_GUI, "Start the IGVtools GUI");
        commandList.put(CMD_HELP, "Usage: igvtools help [command]\n Displays a list of commands."
                + " If [command] is provided, displays help on that command");
    }

    static String[] commandDocs = new String[]{
            "version print the version number",
            "sort    sort an alignment file by start position. ",
            "index   index an alignment file",
            "toTDF    convert an input file (cn, gct, wig) to tiled data format (tdf)",
            "count   compute coverage density for an alignment file",
            "formatexp  center, scale, and log2 normalize an expression file",
            "gui      Start the gui",
            "help     display this help message, or help on a specific command",
            "See http://www.broadinstitute.org/software/igv/igvtools_commandline for more detailed help"
    };
    public static final int MAX_RECORDS_IN_RAM = 500000;
    public static final int MAX_ZOOM = 7;
    public static final int WINDOW_SIZE = 25;
    public static final int EXT_FACTOR = 0;
    public static final Object PROBE_FILE = null;
    public static final int LINEAR_BIN_SIZE = 16000;
    public static final int INTERVAL_SIZE = 1000;
    public static final int LINEAR_INDEX = 1;
    public static final int INTERVAL_INDEX = 2;

    //Options common to all
    private static CmdLineParser.Option windowFunctions = null;
    private static CmdLineParser.Option tmpDirOption = null;
    private static CmdLineParser.Option maxZoomOption = null;
    private static CmdLineParser.Option typeOption = null;

    // options for sort
    private static CmdLineParser.Option maxRecordsOption = null;

    // options for gct files
    private static CmdLineParser.Option probeFileOption = null;

    // options for coverage
    private static CmdLineParser.Option windowSizeOption = null;
    private static CmdLineParser.Option extFactorOption = null;


    private static CmdLineParser.Option separateBasesOption = null;
    private static CmdLineParser.Option strandOption = null;
    private static CmdLineParser.Option queryStringOpt = null;
    private static CmdLineParser.Option minMapQualityOpt = null;
    private static CmdLineParser.Option includeDupsOpt = null;

    // options for index
    private static CmdLineParser.Option indexTypeOption = null;
    private static CmdLineParser.Option binSizeOption = null;
    private static CmdLineParser.Option outputDirOption = null;

    // Trackline
    private static CmdLineParser.Option colorOption = null;


    /**
     * The general usage string
     */
    static String usageString() {
        StringBuffer buf = new StringBuffer();

        buf.append("\nProgram: igvtools\n\n");
        buf.append("Usage: igvtools [command] [options] [input file/dir] [other arguments]\n\n");
        buf.append("Command:");
        for (String c : commandDocs) {
            buf.append(" " + c + "\n\t");
        }
        return buf.toString();
    }

    /**
     * Help on a specific command
     *
     * @param command
     * @return
     */
    static String usageString(String command) {
        command = command.toLowerCase();
        StringBuffer buf = new StringBuffer();
        if (!commandList.containsKey(command)) {
            buf.append("Command " + command + " not recognized\n");
            buf.append(usageString());
            return buf.toString();
        } else {
            buf.append(commandList.get(command));
            return buf.toString();
        }
    }

    private static void initLogger() {
        RollingFileAppender appender = new RollingFileAppender();
        PatternLayout layout = new PatternLayout();
        layout.setConversionPattern("%p [%d{ISO8601}] [%F:%L]  %m%n");
        appender.setName("R");
        appender.setFile("igv.log");
        appender.setThreshold(Level.ALL);
        appender.setMaxFileSize("10KB");
        appender.setMaxBackupIndex(1);
        appender.setAppend(true);
        appender.activateOptions();
        appender.setLayout(layout);
        Logger.getRootLogger().addAppender(appender);
    }


    /**
     * Main method.  Used to run igvtools from the command line
     *
     * @param argv
     * @throws IOException
     * @throws PreprocessingException
     */
    public static void main(String[] argv) throws IOException, PreprocessingException {

        initLogger();
        Globals.setHeadless(true);

        (new IgvTools()).run(argv);

        System.out.println("Done");
        System.exit(1);
    }

    void run(String[] argv) {
        String command = argv[0].toLowerCase();

        if (command.equals(CMD_HELP)) {
            if (argv.length > 1) {
                System.out.println(usageString(argv[1]));
            } else {
                System.out.println(usageString());
            }
        }

        if (command.equals(CMD_GUI)) {
            launchGUI();
            Runtime.getRuntime().halt(0);
        }

        // Do "version" now, its the only command with no arguments
        if (command.equals(CMD_VERSION)) {
            System.out.println("Version " + version);
            return;
        }


        CmdLineParser parser = initParser(command);

        // Parse optional arguments (switches, etc)
        try {
            parser.parse(argv);
        } catch (CmdLineParser.OptionException e) {
            System.err.println(e.getMessage());
            System.out.println(usageString(command));
            return;
        }

        String tmpDirName = null;
        if (tmpDirOption != null) {
            tmpDirName = (String) parser.getOptionValue(tmpDirOption, null);
        }
        int maxRecords = MAX_RECORDS_IN_RAM;
        if (maxRecordsOption != null) {
            maxRecords = (Integer) parser.getOptionValue(maxRecordsOption, MAX_RECORDS_IN_RAM);
        }
        String[] nonOptionArgs = parser.getRemainingArgs();

        try {
            String basic_syntax = "Error in syntax. Expected: " + command + " [options] inputfile outputfile";

            // All remaining commands require an input file, and most need the file extension.  Do that here.
            validateArgsLength(nonOptionArgs, 2, "Error: No input file provided");
            String ifile = nonOptionArgs[1];

            boolean isList = ifile.indexOf(",") > 0;
            if (!isList && !FileUtils.resourceExists(ifile)) {
                throw new PreprocessingException("File not found: " + ifile);
            }

            String typeString = null;
            if (typeOption != null) {
                typeString = (String) parser.getOptionValue(typeOption);
            }
            if (typeString == null || typeString.length() == 0) {
                typeString = Preprocessor.getExtension(ifile).toLowerCase();
            } else {
                typeString = typeString.toLowerCase();
            }


            if (command.equals(CMD_COUNT) || command.equals(CMD_TILE) || command.equals(CMD_TOTDF)) {
                // Parse out options common to both count and tile
                validateArgsLength(nonOptionArgs, 4, commandList.get(command));
                int maxZoomValue = (Integer) parser.getOptionValue(maxZoomOption, MAX_ZOOM);
                String ofile = nonOptionArgs[2];
                String genomeId = nonOptionArgs[3];

                boolean isGCT = typeString.endsWith("gct") || typeString.equals("mage-tab");
                String wfsString = (String) parser.getOptionValue(windowFunctions);
                Collection<WindowFunction> wfList = parseWFS(wfsString, isGCT);


                if (command.equals(CMD_COUNT)) {

                    String trackLine = null;
                    String color = (String) parser.getOptionValue(colorOption);

                    if (color != null) {
                        trackLine = "track color=\"" + color + "\"";
                    }

                    int extFactorValue = (Integer) parser.getOptionValue(extFactorOption, EXT_FACTOR);

                    int countFlags = parseCountFlags(parser);
                    String queryString = (String) parser.getOptionValue(queryStringOpt);
                    int minMapQuality = (Integer) parser.getOptionValue(minMapQualityOpt, 0);

                    int windowSizeValue = (Integer) parser.getOptionValue(windowSizeOption, WINDOW_SIZE);
                    doCount(ifile, ofile, genomeId, maxZoomValue, wfList, windowSizeValue, extFactorValue,
                            trackLine, queryString, minMapQuality, countFlags);
                } else {
                    String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                    toTDF(typeString, ifile, ofile, probeFile, genomeId, maxZoomValue, wfList, tmpDirName, maxRecords);
                }

            } else if (command.equals(CMD_SORT)) {
                validateArgsLength(nonOptionArgs, 3, commandList.get(command));
                String ofile = nonOptionArgs[2];
                doSort(ifile, ofile, tmpDirName, maxRecords);
            } else if (command.equals(CMD_INDEX)) {
                int indexType = (Integer) parser.getOptionValue(indexTypeOption, LINEAR_INDEX);
                int defaultBinSize = indexType == LINEAR_INDEX ? LINEAR_BIN_SIZE : INTERVAL_SIZE;
                int binSize = (Integer) parser.getOptionValue(binSizeOption, defaultBinSize);
                String outputDir = (String) parser.getOptionValue(outputDirOption, null);
                doIndex(ifile, outputDir, indexType, binSize);
            } else if (command.equals(CMD_FORMATEXP)) {
                validateArgsLength(nonOptionArgs, 3, commandList.get(command));
                File inputFile = new File(nonOptionArgs[1]);
                File outputFile = new File(nonOptionArgs[2]);
                (new ExpressionFormatter()).convert(inputFile, outputFile);
            } else if (command.equals("wibtowig")) {
                validateArgsLength(nonOptionArgs, 4, "Error in syntax. Expected: " + command + " [options] txtfile wibfile wigfile");
                File txtFile = new File(nonOptionArgs[1]);
                File wibFile = new File(nonOptionArgs[2]);
                File wigFile = new File(nonOptionArgs[3]);
                String trackLine = nonOptionArgs.length > 4 ? nonOptionArgs[4] : null;
                doWIBtoWIG(txtFile, wibFile, wigFile, trackLine);
            } else if (command.equals("splitgff")) {
                validateArgsLength(nonOptionArgs, 3, "Error in syntax. Expected: " + command + " [options] inputfile outputdir");
                String outputDirectory = nonOptionArgs[2];
                GFFParser.splitFileByType(ifile, outputDirectory);
            } else if (command.toLowerCase().equals("gcttoigv")) {
                validateArgsLength(nonOptionArgs, 4, basic_syntax + " genomeId");
                String ofile = nonOptionArgs[2];
                // Output files must have .igv extension
                if (!ofile.endsWith(".igv")) {
                    ofile = ofile + ".igv";
                }
                String genomeId = nonOptionArgs[3];
                Genome genome = loadGenome(genomeId, true);
                if (genome == null) {
                    throw new PreprocessingException("Genome could not be loaded: " + genomeId);
                }
                String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                doGCTtoIGV(typeString, ifile, new File(ofile), probeFile, maxRecords, tmpDirName, genome);
            } else if (command.equals("formatexp")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                File inputFile = new File(nonOptionArgs[1]);
                File outputFile = new File(nonOptionArgs[2]);
                (new ExpressionFormatter()).convert(inputFile, outputFile);
            } else if (command.toLowerCase().equals("tdftobedgraph")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                String ofile = nonOptionArgs[2];
                TDFUtils.tdfToBedgraph(ifile, ofile);
            } else if (command.equals("wigtobed")) {
                validateArgsLength(nonOptionArgs, 2, "Error in syntax. Expected: " + command + " [options] inputfile");
                String inputFile = nonOptionArgs[1];
                float hetThreshold = 0.17f;
                if (nonOptionArgs.length > 2) {
                    hetThreshold = Float.parseFloat(nonOptionArgs[2]);
                }
                float homThreshold = 0.55f;
                if (nonOptionArgs.length > 3) {
                    homThreshold = Float.parseFloat(nonOptionArgs[3]);
                }
                WigToBed.run(inputFile, hetThreshold, homThreshold);
            } else if (command.equals("vcftobed")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                String inputFile = nonOptionArgs[1];
                String outputFile = nonOptionArgs[2];
                VCFtoBed.convert(inputFile, outputFile);
            } else if (command.equals("lanecounter")) {
                validateArgsLength(nonOptionArgs, 3, "This command is not documented");
                Genome genome = loadGenome(nonOptionArgs[1], false);
                String bamFileList = nonOptionArgs[2];
                String queryInterval = nonOptionArgs[3];
                LaneCounter.run(genome, bamFileList, queryInterval);
            } else if (command.equals("sumwigs")) {
                sumWigs(nonOptionArgs[1], nonOptionArgs[2]);
            } else if (command.equals("densitiestobedgraph")) {
                validateArgsLength(nonOptionArgs, 3, "Error in syntax. Expected: " + command + " [options] inputdir outputdir");
                File inputDir = new File(nonOptionArgs[1]);
                File outputDir = new File(nonOptionArgs[2]);
                if (inputDir.isDirectory() && outputDir.isDirectory()) {
                    DensitiesToBedGraph.convert(inputDir, outputDir);
                } else if (inputDir.isFile() && outputDir.isFile()) {
                    DensitiesToBedGraph.convert(inputDir, outputDir);
                }

            } else {
                throw new PreprocessingException("Unknown command: " + argv[EXT_FACTOR]);
            }
        } catch (PreprocessingException e) {
            System.err.println(e.getMessage());
            System.out.println(usageString());
        } catch (IOException e) {
            throw new PreprocessingException("Unexpected IO error: ", e);
        }
    }


    private CmdLineParser initParser(String command) {
        command = command.toLowerCase();
        CmdLineParser parser = new CmdLineParser();
        if (command.equals(CMD_SORT) || command.equals(CMD_TOTDF) || command.equals(CMD_TILE)) {
            maxRecordsOption = parser.addIntegerOption('m', "maxRecords");
            tmpDirOption = parser.addStringOption('t', "tmpDir");
        }

        if (command.equals(CMD_COUNT) || command.equals(CMD_TOTDF) || command.equals(CMD_TILE)) {

            // general options
            windowFunctions = parser.addStringOption('f', "windowFunctions");
            maxZoomOption = parser.addIntegerOption('z', "maxZoom");

            // extended options for coverage
            if (command.equals(CMD_COUNT)) {

                extFactorOption = parser.addIntegerOption('e', "extFactor");
                windowSizeOption = parser.addIntegerOption('w', "windowSize");

                separateBasesOption = parser.addBooleanOption("bases");
                strandOption = parser.addStringOption("strand");
                queryStringOpt = parser.addStringOption("query");
                minMapQualityOpt = parser.addStringOption("minMapQuality");
                includeDupsOpt = parser.addStringOption("includeDuplicates");

                // Trackline
                colorOption = parser.addStringOption("color");
            } else {
                // options for gct files
                probeFileOption = parser.addStringOption('p', "probeFile");
                typeOption = parser.addStringOption("fileType");
            }
        }

        if (command.equals(CMD_INDEX)) {
            indexTypeOption = parser.addIntegerOption("indexType");
            binSizeOption = parser.addIntegerOption("binSize");
            outputDirOption = parser.addStringOption("outputDir");
        }

        return parser;
    }

    private int parseCountFlags(CmdLineParser parser) {

        int countFlags = 0;
        countFlags += (Boolean) parser.getOptionValue(separateBasesOption, false) ? CoverageCounter.BASES : 0;
        countFlags += (Boolean) parser.getOptionValue(includeDupsOpt, false) ? CoverageCounter.INCLUDE_DUPS : 0;
        String strandopt = (String) parser.getOptionValue(strandOption, "");
        if (strandopt.equals("read")) {
            countFlags += CoverageCounter.STRAND_SEPARATE;
        } else if (strandopt.equals("first")) {
            countFlags += CoverageCounter.FIRST_IN_PAIR;
        }
        return countFlags;
    }

    void run_old(String[] argv) {

        CmdLineParser parser = new CmdLineParser();
        CmdLineParser.Option helpOption = parser.addBooleanOption('h', "help");
        CmdLineParser.Option guiOption = parser.addBooleanOption('g', "gui");

        // general options
        CmdLineParser.Option windowFunctions = parser.addStringOption('f', "windowFunctions");
        CmdLineParser.Option tmpDirOption = parser.addStringOption('t', "tmpDir");
        CmdLineParser.Option maxZoomOption = parser.addIntegerOption('z', "maxZoom");

        // options for sort
        CmdLineParser.Option maxRecordsOption = parser.addIntegerOption('m', "maxRecords");

        // options for gct files
        CmdLineParser.Option probeFileOption = parser.addStringOption('p', "probeFile");

        // options for coverage
        CmdLineParser.Option windowSizeOption = parser.addIntegerOption('w', "windowSize");
        CmdLineParser.Option extFactorOption = parser.addIntegerOption('e', "extFactor");

        // options for index
        CmdLineParser.Option indexTypeOption = parser.addIntegerOption('i', "indexType");
        CmdLineParser.Option binSizeOption = parser.addIntegerOption('b', "binSize");
        CmdLineParser.Option outputDirOption = parser.addStringOption('o', "outputDir");

        // extended options for coverage
        //  q, t, .....
        CmdLineParser.Option coverageOptions = parser.addStringOption('a', "coverageOptions");

        // Trackline
        CmdLineParser.Option colorOption = parser.addStringOption('c', "color");


        // Parse optional arguments (switches, etc)
        try {
            parser.parse(argv);
        } catch (CmdLineParser.OptionException e) {
            System.err.println(e.getMessage());
            return;
        }

        boolean help = ((Boolean) parser.getOptionValue(helpOption, false));
        if (help) {
            System.out.println(usageString());
            return;
        }


        boolean gui = ((Boolean) parser.getOptionValue(guiOption, false));
        if (gui) {
            launchGUI();
            Runtime.getRuntime().halt(0);
        }


        String tmpDirName = (String) parser.getOptionValue(tmpDirOption);


        String[] nonOptionArgs = parser.getRemainingArgs();

        try {
            validateArgsLength(argv, 1, "Error: Must have at least 1 argument");
            // The command
            String command = nonOptionArgs[EXT_FACTOR].toLowerCase();
            String basic_syntax = "Error in syntax. Expected: " + command + " [options] inputfile outputfile";

            // Do "version" now, its the only command with no arguments
            if (command.equals("version")) {
                System.out.println("Version " + version);
                return;
            }

            // All remaining commands require an input file, and most need the file extension.  Do that here.
            validateArgsLength(nonOptionArgs, 2, "Error: No input file provided");
            String ifile = nonOptionArgs[1];
            String typeString = (String) parser.getOptionValue(typeOption);
            if (typeString == null) {
                typeString = Preprocessor.getExtension(ifile).toLowerCase();
            } else {
                typeString = typeString.toLowerCase();
            }

            boolean isList = ifile.indexOf(",") > 0;
            if (!isList && !FileUtils.resourceExists(ifile)) {
                throw new PreprocessingException("File not found: " + ifile);
            }

            int maxRecords = (Integer) parser.getOptionValue(maxRecordsOption, MAX_RECORDS_IN_RAM);

            if (command.equals("count") || command.equals("totdf") || command.equals("tile")) {

                // Parse out options common to both count and tile
                validateArgsLength(nonOptionArgs, 4, basic_syntax + " genomeId");
                int windowSizeValue = (Integer) parser.getOptionValue(windowSizeOption, WINDOW_SIZE);
                int maxZoomValue = (Integer) parser.getOptionValue(maxZoomOption, MAX_ZOOM);
                String ofile = nonOptionArgs[2];
                String genomeId = nonOptionArgs[3];
                boolean isGCT = typeString.endsWith("gct") || typeString.equals("mage-tab");
                String wfsString = (String) parser.getOptionValue(windowFunctions);
                Collection<WindowFunction> wfList = parseWFS(wfsString, isGCT);


                String coverageOpt = (String) parser.getOptionValue(coverageOptions);

                String trackLine = null;
                String color = (String) parser.getOptionValue(colorOption);
                if (color != null) {
                    trackLine = "track color=\"" + color + "\"";
                }

                if (command.equals("count")) {
                    int extFactorValue = (Integer) parser.getOptionValue(extFactorOption, EXT_FACTOR);

                    //doCount(ifile, ofile, genomeId, maxZoomValue, wfList, windowSizeValue, extFactorValue,
                    //coverageOpt, trackLine);
                } else {
                    String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                    toTDF(typeString, ifile, ofile, probeFile, genomeId, maxZoomValue, wfList, tmpDirName, maxRecords);
                }
            } else if (command.equals("sort")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);

                String ofile = nonOptionArgs[2];
                doSort(ifile, ofile, tmpDirName, maxRecords);
            } else if (command.equals("index")) {
                int indexType = (Integer) parser.getOptionValue(indexTypeOption, LINEAR_INDEX);
                int defaultBinSize = indexType == LINEAR_INDEX ? LINEAR_BIN_SIZE : INTERVAL_SIZE;
                int binSize = (Integer) parser.getOptionValue(binSizeOption, defaultBinSize);
                String outputDir = (String) parser.getOptionValue(outputDirOption, null);
                doIndex(ifile, outputDir, indexType, binSize);
            } else if (command.equals("wibtowig")) {
                validateArgsLength(nonOptionArgs, 4, "Error in syntax. Expected: " + command + " [options] txtfile wibfile wigfile");
                File txtFile = new File(nonOptionArgs[1]);
                File wibFile = new File(nonOptionArgs[2]);
                File wigFile = new File(nonOptionArgs[3]);
                String trackLine = nonOptionArgs.length > 4 ? nonOptionArgs[4] : null;
                doWIBtoWIG(txtFile, wibFile, wigFile, trackLine);
            } else if (command.equals("splitgff")) {
                validateArgsLength(nonOptionArgs, 3, "Error in syntax. Expected: " + command + " [options] inputfile outputdir");
                String outputDirectory = nonOptionArgs[2];
                GFFParser.splitFileByType(ifile, outputDirectory);
            } else if (command.toLowerCase().equals("gcttoigv")) {
                validateArgsLength(nonOptionArgs, 4, basic_syntax + " genomeId");
                String ofile = nonOptionArgs[2];
                // Output files must have .igv extension
                if (!ofile.endsWith(".igv")) {
                    ofile = ofile + ".igv";
                }
                String genomeId = nonOptionArgs[3];
                Genome genome = loadGenome(genomeId, true);
                if (genome == null) {
                    throw new PreprocessingException("Genome could not be loaded: " + genomeId);
                }
                String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                doGCTtoIGV(typeString, ifile, new File(ofile), probeFile, maxRecords, tmpDirName, genome);
            } else if (command.equals("formatexp")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                File inputFile = new File(nonOptionArgs[1]);
                File outputFile = new File(nonOptionArgs[2]);
                (new ExpressionFormatter()).convert(inputFile, outputFile);
            } else if (command.toLowerCase().equals("tdftobedgraph")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                String ofile = nonOptionArgs[2];
                TDFUtils.tdfToBedgraph(ifile, ofile);
            } else if (command.equals("wigtobed")) {
                validateArgsLength(nonOptionArgs, 2, "Error in syntax. Expected: " + command + " [options] inputfile");
                String inputFile = nonOptionArgs[1];
                float hetThreshold = 0.17f;
                if (nonOptionArgs.length > 2) {
                    hetThreshold = Float.parseFloat(nonOptionArgs[2]);
                }
                float homThreshold = 0.55f;
                if (nonOptionArgs.length > 3) {
                    homThreshold = Float.parseFloat(nonOptionArgs[3]);
                }
                WigToBed.run(inputFile, hetThreshold, homThreshold);
            } else if (command.equals("vcftobed")) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                String inputFile = nonOptionArgs[1];
                String outputFile = nonOptionArgs[2];
                VCFtoBed.convert(inputFile, outputFile);
            } else if (command.equals("lanecounter")) {
                //TODO
                validateArgsLength(nonOptionArgs, 3, "This command is not documented");
                Genome genome = loadGenome(nonOptionArgs[1], false);
                String bamFileList = nonOptionArgs[2];
                String queryInterval = nonOptionArgs[3];
                LaneCounter.run(genome, bamFileList, queryInterval);
            } else if (command.equals("sumwigs")) {
                sumWigs(nonOptionArgs[1], nonOptionArgs[2]);
            } else if (command.equals("densitiestobedgraph")) {
                validateArgsLength(nonOptionArgs, 3, "Error in syntax. Expected: " + command + " [options] inputdir outputdir");
                File inputDir = new File(nonOptionArgs[1]);
                File outputDir = new File(nonOptionArgs[2]);
                if (inputDir.isDirectory() && outputDir.isDirectory()) {
                    DensitiesToBedGraph.convert(inputDir, outputDir);
                } else if (inputDir.isFile() && outputDir.isFile()) {
                    DensitiesToBedGraph.convert(inputDir, outputDir);
                }

            } else {
                throw new PreprocessingException("Unknown command: " + argv[EXT_FACTOR]);
            }
        } catch (PreprocessingException e) {
            System.err.println(e.getMessage());
        } catch (IOException e) {
            throw new PreprocessingException("Unexpected IO error: ", e);
        }
    }

    private void sumWigs(String inputString, String outputString) throws IOException {

        String[] tokens = inputString.split(",");
        List<File> in = new ArrayList();
        for (String f : tokens) {
            in.add(new File(f));
        }
        File out = new File(outputString);
        WigSummer.sumWigs(in, out);

    }

    private void doGCTtoIGV(String typeString, String ifile, File ofile, String probefile, int maxRecords, String tmpDirName, Genome genome) throws IOException {

        System.out.println("gct -> igv: " + ifile + " -> " + ofile.getAbsolutePath());

        File tmpDir = null;
        if (tmpDirName != null && tmpDirName.trim().length() > 0) {
            tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                System.err.println("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
                throw new PreprocessingException("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
            }
        }

        ResourceLocator locator = new ResourceLocator(ifile);
        locator.setType(typeString);
        GCTtoIGVConverter.convert(locator, ofile, probefile, maxRecords, tmpDir, genome);

    }

    public void toTDF(String typeString, String ifile, String ofile, String probeFile, String genomeId, int maxZoomValue,
                      Collection<WindowFunction> windowFunctions, String tmpDirName, int maxRecords)
            throws IOException, PreprocessingException {

        if (!ifile.endsWith(".affective.csv")) validateIsTilable(typeString);

        System.out.println("Tile.  File = " + ifile);
        System.out.println("Max zoom = " + maxZoomValue);
        if (probeFile != null && probeFile.trim().length() > 0) {
            System.out.println("Probe file = " + probeFile);
        }
        System.out.print("Window functions: ");
        for (WindowFunction wf : windowFunctions) {
            System.out.print(wf.toString() + " ");
        }
        System.out.println();

        boolean isGCT = isGCT(typeString);
        Genome genome = loadGenome(genomeId, isGCT);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }
        File tmp = new File(ifile);


        int nLines = tmp.isDirectory() ? ParsingUtils.estimateLineCount(tmp) : ParsingUtils.estimateLineCount(ifile);

        // Convert  gct files to igv format first
        File deleteme = null;
        if (isGCT(typeString)) {
            File tmpDir = null;
            if (tmpDirName != null && tmpDirName.length() > 0) {
                tmpDir = new File(tmpDirName);
                if (!tmpDir.exists() || !tmpDir.isDirectory()) {
                    throw new PreprocessingException("Specified tmp directory does not exist or is not directory: " + tmpDirName);
                }
            } else {
                tmpDir = new File(System.getProperty("java.io.tmpdir"), System.getProperty("user.name"));
            }
            if (!tmpDir.exists()) {
                tmpDir.mkdir();
            }

            String baseName = (new File(ifile)).getName();
            File igvFile = new File(tmpDir, baseName + ".igv");
            igvFile.deleteOnExit();
            doGCTtoIGV(typeString, ifile, igvFile, probeFile, maxRecords, tmpDirName, genome);

            tmp = igvFile;
            deleteme = igvFile;

        }


        File outputFile = new File(ofile);
        try {
            Preprocessor p = new Preprocessor(outputFile, genome, windowFunctions, nLines, null);
            if (tmp.isDirectory()) {
                File[] files = tmp.listFiles();
                Arrays.sort(files, new Comparator<File>() {
                    public int compare(File file, File file1) {
                        return file.getName().compareTo(file1.getName());
                    }
                });
                for (File f : files) {
                    p.preprocess(f, maxZoomValue);
                }
            } else {
                p.preprocess(tmp, maxZoomValue);
            }
            p.finish();
        } catch (IOException e) {
            e.printStackTrace();
            // Delete output file as its probably corrupt
            if (outputFile.exists()) {
                outputFile.delete();
            }
        } finally {
            if (deleteme != null && deleteme.exists()) {
                deleteme.delete();
            }
        }

        System.out.flush();

    }

    private boolean isGCT(String tmp) {
        if (tmp.endsWith(".txt")) tmp = tmp.substring(0, tmp.length() - 4);
        if (tmp.endsWith(".gz")) tmp = tmp.substring(0, tmp.length() - 3);
        return tmp.endsWith(".gct") || tmp.endsWith(".tab") || tmp.equals("mage-tab");
    }

    /**
     * Compute coverage or density of an alignment or feature file.
     *
     * @param ifile           Alignement or feature file
     * @param ofile           Output file
     * @param genomeId        Genome id (e.g. hg18) or full path to a .genome file (e.g. /xchip/igv/scer2.genome)
     * @param maxZoomValue    Maximum zoom level to precompute.  Default value is 7
     * @param windowFunctions
     * @param windowSizeValue
     * @param extFactorValue
     * @throws IOException
     */
    public void doCount(String ifile, String ofile, String genomeId, int maxZoomValue,
                        Collection<WindowFunction> windowFunctions, int windowSizeValue, int extFactorValue,
                        String trackLine, String queryString, int minMapQuality, int countFlags) throws IOException {


        System.out.println("Computing coverage.  File = " + ifile);
        System.out.println("Max zoom = " + maxZoomValue);
        System.out.println("Window size = " + windowSizeValue);
        System.out.print("Window functions: ");
        for (WindowFunction wf : windowFunctions) {
            System.out.print(wf.toString() + " ");
        }
        System.out.println();
        System.out.println("Ext factor = " + extFactorValue);


        Genome genome = loadGenome(genomeId, false);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }

        // Multiple files allowed for count command (a tdf and a wig)
        File tdfFile = null;
        File wigFile = null;
        String[] files = ofile.split(",");
        if (files[0].endsWith("wig")) {
            wigFile = new File(files[0]);
        } else {
            tdfFile = new File(files[0]);
        }
        if (files.length > 1) {
            if (files[1].endsWith("wig")) {
                wigFile = new File(files[1]);
            } else if (files[1].endsWith("tdf")) {
                tdfFile = new File(files[1]);
            }
        }

        if (tdfFile != null && !tdfFile.getName().endsWith(".tdf")) {
            tdfFile = new File(tdfFile.getAbsolutePath() + ".tdf");
        }

        Preprocessor p = new Preprocessor(tdfFile, genome, windowFunctions, -1, null);
        //p.count(ifile, windowSizeValue, extFactorValue, maxZoomValue, wigFile, coverageOpt, trackLine);
        p.count(ifile, windowSizeValue, extFactorValue, maxZoomValue, wigFile, trackLine,
                queryString, minMapQuality, countFlags);

        p.finish();

        System.out.flush();
    }


    public void doWIBtoWIG(File txtFile, File wibFile, File wigFile, String trackLine) {
        UCSCUtils.convertWIBFile(txtFile, wibFile, wigFile, trackLine);
    }

    /**
     * Create an index for an alignment or feature file
     *
     * @param ifile
     * @throws IOException
     */

    public void doIndex(String ifile, int indexType, int binSize) throws IOException {
        doIndex(ifile, null, indexType, binSize);
    }

    public void doIndex(String ifile, String outputFileName, int indexType, int binSize) throws IOException {
        if (ifile.endsWith(".gz")) {
            System.out.println("Cannot index a gzipped file");
            throw new PreprocessingException("Cannot index a gzipped file");
        }

        if (outputFileName == null) {
            outputFileName = ifile;
        }

        FeatureCodec codec = CodecFactory.getCodec(ifile);
        if (codec != null) {

            if (!outputFileName.endsWith(".idx")) {
                outputFileName = outputFileName + ".idx";
            }
            try {
                createTribbleIndex(ifile, new File(outputFileName), indexType, binSize, codec);
            } catch (TribbleException.MalformedFeatureFile e) {
                StringBuffer buf = new StringBuffer();
                buf.append("<html>Files must be sorted by start position prior to indexing.<br>");
                buf.append(e.getMessage());
                buf.append("<br><br>Note: igvtools can be used to sort the file, select \"File > Run igvtools...\".");
                MessageUtils.showMessage(buf.toString());
            }


        } else {
            if (!outputFileName.endsWith(".sai")) {
                outputFileName = outputFileName + ".sai";
            }
            AlignmentIndexer indexer = AlignmentIndexer.getInstance(new File(ifile), null, null);
            File outputFile = new File(outputFileName);
            try {
                indexer.createSamIndex(outputFile);
            } catch (Exception e) {
                e.printStackTrace();
                // Delete output file as it is probably corrupt
                if (outputFile.exists()) {
                    outputFile.delete();
                }
            }
        }
        System.out.flush();

    }

    /**
     * Create a tribble style index.
     *
     * @param ifile
     * @param outputFile
     * @param indexType
     * @param binSize
     * @throws IOException
     */
    private void createTribbleIndex(String ifile, File outputFile, int indexType, int binSize, FeatureCodec codec) throws IOException {
        File inputFile = new File(ifile);
        Index idx = null;
        if (indexType == LINEAR_INDEX) {
            idx = IndexFactory.createLinearIndex(inputFile, codec, binSize);
        } else {
            idx = IndexFactory.createIntervalIndex(inputFile, codec, binSize);
        }
        if (idx != null) {
            String idxFile;
            if (outputFile != null) {
                idxFile = outputFile.getAbsolutePath();
            } else {
                idxFile = ifile = ".idx";
            }
            LittleEndianOutputStream stream = null;
            try {
                stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
                idx.write(stream);
            } catch (Exception e) {
                e.printStackTrace();
                // Delete output file as its probably corrupt
                File tmp = new File(idxFile);
                if (tmp.exists()) {
                    tmp.delete();
                }
            } finally {
                if (stream != null) {
                    stream.close();
                }
            }
        }
    }


    public void doSort(String ifile, String ofile, String tmpDirName, int maxRecords) {

        System.out.println("Sorting " + ifile + "  -> " + ofile);
        File inputFile = new File(ifile);
        File outputFile = new File(ofile);
        Sorter sorter = Sorter.getSorter(inputFile, outputFile);
        if (tmpDirName != null && tmpDirName.trim().length() > 0) {
            File tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                System.err.println("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
                throw new PreprocessingException("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
            }
            sorter.setTmpDir(tmpDir);
        }

        sorter.setMaxRecords(maxRecords);
        try {
            sorter.run();
        } catch (Exception e) {
            e.printStackTrace();
            // Delete output file as its probably corrupt
            if (outputFile.exists()) {
                outputFile.delete();
            }
        }
        System.out.println("Done");
        System.out.flush();
    }


    private void validateArgsLength(String[] nonOptionArgs, int len, String failMessage) throws PreprocessingException {
        if (nonOptionArgs.length < len) {
            throw new PreprocessingException(failMessage + "\n" + usageString());
        }
    }


    public static Genome loadGenome(String genomeFileOrID, boolean loadGenes) throws IOException {

        String rootDir = FileUtils.getInstallDirectory();

        final GenomeManager genomeManager = Globals.isHeadless() ? new GenomeManager() :
                IGV.getInstance().getGenomeManager();
        Genome genome = genomeManager.getCurrentGenome();
        if (genome != null && genome.getId().equals(genomeFileOrID)) {
            return genome;
        }

        File genomeFile = new File(genomeFileOrID);
        if (!genomeFile.exists()) {
            genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID + ".genome");

        }
        if (!genomeFile.exists()) {
            genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID);

        }
        if (!genomeFile.exists()) {
            genomeFile = new File(genomeFileOrID);
        }
        if (!genomeFile.exists()) {
            throw new PreprocessingException("Genome definition file not found for: " + genomeFileOrID);
        }

        genome = genomeManager.loadGenome(genomeFile.getAbsolutePath(), null);
        if (genome == null) {
            throw new PreprocessingException("Error loading: " + genomeFileOrID);
        }

        // If this is a .genome file optionally file load genes
        if (loadGenes && genomeFile.getAbsolutePath().endsWith(".genome")) {
            GenomeDescriptor descriptor = genomeManager.parseGenomeArchiveFile(genomeFile);
            String geneFileName = descriptor.getGeneFileName();
            if (geneFileName != null && geneFileName.trim().length() > 0) {
                FeatureParser parser = AbstractFeatureParser.getInstanceFor(geneFileName, genome);
                InputStream is = descriptor.getGeneStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));
                //Right now the parser adds these to the FeatureDB map
                //May want to move that someplace else
                List<Feature> features = parser.loadFeatures(reader);
                is.close();
            }
        }


        return genome;
    }


    /**
     * Test if the file type can be "tiled".
     */
    private static void validateIsTilable(String typeString) {

        boolean affective = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.AFFECTIVE_ENABLE);
        if (!(typeString.equals(".cn") ||
                typeString.equals(".igv") ||
                typeString.equals(".wig") ||
                // ifile.toLowerCase().endsWith("cpg.txt") ||
                typeString.equals(".ewig") ||
                typeString.equals(".cn") ||
                typeString.equals(".snp") ||
                typeString.equals(".xcn") ||
                typeString.equals(".gct") ||
                typeString.equals("mage-tab") ||
                typeString.equals(".bedgraph") ||
                Preprocessor.isAlignmentFile(typeString) ||
                affective)) {
            throw new PreprocessingException("Tile command not supported for files of type: " + typeString);
        }
    }


    /**
     * Parse the window functions line.   The default for most files is a single "mean",  however gct files include
     * min and max as well.
     *
     * @param string comma delimited string of window functions, e.g. min, p10, max
     * @return colleciton of WindowFunctions objects
     */
    private static Collection<WindowFunction> parseWFS(String string, boolean isGCT) {
        if (string == null || string.length() == EXT_FACTOR) {
            return isGCT ? Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max) :
                    Arrays.asList(WindowFunction.mean);
        } else {
            String[] tokens = string.split(",");
            List<WindowFunction> funcs = new ArrayList(tokens.length);
            for (int i = EXT_FACTOR; i < tokens.length; i++) {
                String wf = tokens[i];
                if (wf.startsWith("p")) {
                    wf = wf.replaceFirst("p", "percentile");
                }
                try {
                    funcs.add(WindowFunction.valueOf(wf));
                } catch (Exception e) {
                    System.err.println("Unrecognized window function: " + tokens[i]);
                }
            }
            return funcs;
        }
    }


    private static void launchGUI() {
        IgvToolsGui.main(null);
    }

}

