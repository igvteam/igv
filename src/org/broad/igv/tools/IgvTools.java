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
package org.broad.igv.tools;


import jargs.gnu.CmdLineParser;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.*;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.genome.FastaUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDescriptor;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.sam.reader.AlignmentIndexer;
import org.broad.igv.tdf.TDFUtils;
import org.broad.igv.tools.converters.BamToBed;
import org.broad.igv.tools.converters.ExpressionFormatter;
import org.broad.igv.tools.converters.GCTtoIGVConverter;
import org.broad.igv.tools.converters.WigToBed;
import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.ReadmeParser;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.converters.DensitiesToBedGraph;
import org.broad.igv.variant.util.VCFtoBed;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;

/**
 * Command accessories for IGV.
 *
 * @author jrobinso
 */
public class IgvTools {

    static private Logger log = Logger.getLogger(IgvTools.class);

    static final String CMD_TILE = "tile";
    static final String CMD_TOTDF = "totdf";
    static final String CMD_COUNT = "count";
    static final String CMD_SORT = "sort";
    static final String CMD_INDEX = "index";
    static final String CMD_FORMATEXP = "formatexp";
    static final String CMD_VERSION = "version";
    static final String CMD_GUI = "gui";
    static final String CMD_HELP = "help";
    static final String CMD_BAMTOBED = "bamtobed";
    static final String CMD_TDFTOBEDGRAPH = "tdftobedgraph";

    /**
     * Stream for writing messages to the user, which we
     * DO NOT want to permanently log. Anything we want to
     * permanently log should use {@link #log}
     */
    static PrintStream userMessageWriter = System.out;
    private static final String CONSOLE_APPENDER_NAME = "console";

    /**
     * String used in place of output file path to indicate output
     * should be to stdout
     */
    static final String STDOUT_FILE_STR = "stdout";

    public static String getVersionString() {
        return Globals.applicationString();
    }

    //TODO extract this from readme
    static String[] commandDocs = new String[]{
            "version print the version number",
            "sort    sort an alignment file by start position. ",
            "index   index an alignment file",
            "toTDF    convert an input file (cn, gct, wig) to tiled data format (tdf)",
            "count   compute coverage density for an alignment file",
            "formatexp  center, scale, and log2 normalize an expression file",
            "gui      Start the gui",
            "help <command>     display this help message, or help on a specific command",
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
    private static CmdLineParser.Option preExtFactorOption = null;
    private static CmdLineParser.Option postExtFactorOption = null;


    private static CmdLineParser.Option separateBasesOption = null;
    private static CmdLineParser.Option strandOption = null;
    private static CmdLineParser.Option queryStringOpt = null;
    private static CmdLineParser.Option minMapQualityOpt = null;
    private static CmdLineParser.Option includeDupsOpt = null;
    private static CmdLineParser.Option pairedCoverageOpt = null;

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

        buf.append("\nProgram: igvtools. " + getVersionString() + "\n\n");
        buf.append("Usage: igvtools [command] [options] [input file/dir] [other arguments]\n\n");
        buf.append("Command:");
        for (String c : commandDocs) {
            buf.append(" " + c + "\n\t");
        }
        return buf.toString();
    }

    /**
     * Help on a specific command.
     * Parse the README file.
     *
     * @param command
     * @return
     */
    static String usageString(String command) {
        command = command.toLowerCase();
        ReadmeParser parser = new ReadmeParser();
        String help = parser.getDocForCommand(command);
        return help;
    }

    private static void initLogger() {
        if(Logger.getRootLogger().getAppender("R") == null){
            RollingFileAppender fileAppender = new RollingFileAppender();
            PatternLayout fileLayout = new PatternLayout();
            fileLayout.setConversionPattern("%p [%d{ISO8601}] [%F:%L]  %m%n");
            fileAppender.setName("R");
            fileAppender.setFile("igv.log");
            fileAppender.setThreshold(Level.ALL);
            fileAppender.setMaxFileSize("10KB");
            fileAppender.setMaxBackupIndex(1);
            fileAppender.setAppend(true);
            fileAppender.activateOptions();
            fileAppender.setLayout(fileLayout);

            //Logger.getRootLogger().addAppender(fileAppender);
        }

        if(Logger.getRootLogger().getAppender(CONSOLE_APPENDER_NAME) == null){
            PatternLayout consoleLayout = new PatternLayout();
            consoleLayout.setConversionPattern("%m%n");
            ConsoleAppender consoleAppender = new ConsoleAppender();
            consoleAppender.setThreshold(Level.INFO);
            consoleAppender.setName(CONSOLE_APPENDER_NAME);
            consoleAppender.setFollow(true);
            consoleAppender.activateOptions();
            consoleAppender.setLayout(consoleLayout);

            Logger.getRootLogger().addAppender(consoleAppender);
        }
    }


    /**
     * Main method.  Used to run igvtools from the command line
     *
     * @param argv
     * @throws IOException
     * @throws PreprocessingException
     */
    public static void main(String[] argv) {

        try {
            Globals.setHeadless(true);

            (new IgvTools()).run(argv);

            userMessageWriter.println("Done");
            System.exit(0);
        } catch (Exception e) {
            log.error(e.getMessage(), e);
            System.exit(-1);
        }
    }

    void run(String[] argv) {
        initLogger();

        if (argv.length == 0) {
            userMessageWriter.println(usageString());
            userMessageWriter.println("Error: No arguments provided");
            return;
        }

        String command = argv[0].toLowerCase();

        if (command.equals(CMD_HELP)) {
            if (argv.length > 1) {
                userMessageWriter.println(usageString(argv[1]));
            } else {
                userMessageWriter.println(usageString());
            }
            return;
        }

        if (command.equals(CMD_GUI)) {
            launchGUI();
            Runtime.getRuntime().halt(0);
        }

        // Do "version" now, its the only command with no arguments
        if (command.equals(CMD_VERSION)) {
            userMessageWriter.println(getVersionString());
            return;
        }

        CmdLineParser parser = initParser(command);

        // Parse optional arguments (switches, etc)
        try {
            parser.parse(argv);
        } catch (CmdLineParser.OptionException e) {
            userMessageWriter.println(e.getMessage());
            userMessageWriter.println("Enter igvtools help " + command + " for help on this command");
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
            String basic_syntax = "Error in syntax. Enter igvtools help " + command + " for usage instructions.";

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
                validateArgsLength(nonOptionArgs, 4, basic_syntax);
                int maxZoomValue = (Integer) parser.getOptionValue(maxZoomOption, MAX_ZOOM);
                String ofile = nonOptionArgs[2];
                setWriteToStdOout(ofile);

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
                    int preFactorValue = (Integer) parser.getOptionValue(preExtFactorOption, 0);
                    int posFactorValue = (Integer) parser.getOptionValue(postExtFactorOption, 0);

                    int countFlags = parseCountFlags(parser);
                    String queryString = (String) parser.getOptionValue(queryStringOpt);
                    int minMapQuality = (Integer) parser.getOptionValue(minMapQualityOpt, 0);

                    int windowSizeValue = (Integer) parser.getOptionValue(windowSizeOption, WINDOW_SIZE);
                    doCount(ifile, ofile, genomeId, maxZoomValue, wfList, windowSizeValue, extFactorValue,
                            preFactorValue, posFactorValue,
                            trackLine, queryString, minMapQuality, countFlags);
                } else {
                    String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                    toTDF(typeString, ifile, ofile, probeFile, genomeId, maxZoomValue, wfList, tmpDirName, maxRecords);
                }

            } else if (command.equals(CMD_SORT)) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                String ofile = nonOptionArgs[2];
                doSort(ifile, ofile, tmpDirName, maxRecords);
            } else if (command.equals(CMD_INDEX)) {
                int indexType = (Integer) parser.getOptionValue(indexTypeOption, LINEAR_INDEX);
                int defaultBinSize = indexType == LINEAR_INDEX ? LINEAR_BIN_SIZE : INTERVAL_SIZE;
                int binSize = (Integer) parser.getOptionValue(binSizeOption, defaultBinSize);
                String outputDir = (String) parser.getOptionValue(outputDirOption, null);
                doIndex(ifile, typeString, outputDir, indexType, binSize);
            } else if (command.equals(CMD_FORMATEXP)) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
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
            } else if(command.equals("gfftobed")){
                validateArgsLength(nonOptionArgs, 3, "Error in syntax. Expected: " + command + " inputfile outputfile");
                String ofile = nonOptionArgs[2];
                setWriteToStdOout(ofile);
                GFFToBed(ifile, ofile);
            } else if (command.toLowerCase().equals("gcttoigv")) {
                validateArgsLength(nonOptionArgs, 4, basic_syntax + " genomeId");
                String ofile = nonOptionArgs[2];
                // Output files must have .igv extension
                if (!ofile.endsWith(".igv")) {
                    ofile = ofile + ".igv";
                }
                String genomeId = nonOptionArgs[3];
                Genome genome = loadGenome(genomeId);
                if (genome == null) {
                    throw new PreprocessingException("Genome could not be loaded: " + genomeId);
                }
                String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                doGCTtoIGV(typeString, ifile, new File(ofile), probeFile, maxRecords, tmpDirName, genome);
            } else if (command.toLowerCase().equals(CMD_TDFTOBEDGRAPH)) {
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

            } else if (command.equals(CMD_BAMTOBED)) {
                validateArgsLength(nonOptionArgs, 3, basic_syntax);
                String ofile = nonOptionArgs[2];
                Boolean pairOption = (Boolean) parser.getOptionValue(pairedCoverageOpt, false);
                BamToBed.convert(new File(ifile), new File(ofile), pairOption);
            } else if (command.equalsIgnoreCase("genGenomeList")) {
                //Generate a genomes.txt list file based on a directory
                //TODO Probably a better place for this. Users won't generally use it
                File inDir = new File(ifile);
                GenomeManager manager = GenomeManager.getInstance();
                manager.generateGenomeList(inDir, nonOptionArgs[2], nonOptionArgs[3]);
            } else {
                throw new PreprocessingException("Unknown command: " + argv[EXT_FACTOR]);
            }
        } catch (IOException e) {
            throw new PreprocessingException("Unexpected IO error: ", e);
        }
    }

    private void GFFToBed(String ifile, String ofile) throws FileNotFoundException{
        IGVBEDCodec outCodec = new IGVBEDCodec();
        GFFParser parser = new GFFParser();
        GFFCodec codec = null;
        try {
            codec = (GFFCodec) CodecFactory.getCodec(ifile, null);
        } catch (Exception e) {
            throw new IllegalArgumentException("Input file is not recognized as a GFF");
        }
        BufferedReader reader = null;
        PrintStream outStream = System.out;
        if(!ofile.equals(STDOUT_FILE_STR)){
            outStream = new PrintStream(new FileOutputStream(ofile));
        }
        try {
            reader = ParsingUtils.openBufferedReader(ifile);
            List<Feature> features = parser.loadFeatures(reader, null, codec);
            for (Feature feat : features) {
                String encoded = outCodec.encode(feat);
                outStream.print(encoded);
                outStream.print('\n');
            }
        } catch (IOException e) {
            log.error(e.getMessage(), e);
        } finally {
            if (reader != null){
                try {
                    reader.close();
                } catch (IOException e) {
                    log.error(e.getMessage(), e);
                }
            }
            if (outStream != null) {
                outStream.flush();
                outStream.close();
            }
        }
    }

    /**
     * if ofile.equals(STDOUT_FILE_STR), write output to stdout. This also means redirecting log statements
     * to someplace other than stdout, we use stderr
     * @param ofile
     * @return Whether output will be written to stdout
     */
    private boolean setWriteToStdOout(String ofile) {
        //Output will be written to stdout instead of file,
        //need to redirect user messages
        if(ofile.equals(STDOUT_FILE_STR)){

            userMessageWriter = System.err;

            ConsoleAppender appender = (ConsoleAppender) Logger.getRootLogger().getAppender(CONSOLE_APPENDER_NAME);
            appender.setTarget(ConsoleAppender.SYSTEM_ERR);
            appender.activateOptions();

            //See log4j.properties file
            Logger.getRootLogger().removeAppender("stdout");
            return true;
        }
        return false;
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
            if (command.equals(CMD_COUNT) || command.equals(CMD_BAMTOBED)) {

                extFactorOption = parser.addIntegerOption('e', "extFactor");
                preExtFactorOption = parser.addIntegerOption("preExtFactor");
                postExtFactorOption = parser.addIntegerOption("postExtFactor");
                windowSizeOption = parser.addIntegerOption('w', "windowSize");

                separateBasesOption = parser.addBooleanOption("bases");
                strandOption = parser.addStringOption("strands");
                queryStringOpt = parser.addStringOption("query");
                minMapQualityOpt = parser.addIntegerOption("minMapQuality");
                includeDupsOpt = parser.addBooleanOption("includeDuplicates");
                pairedCoverageOpt = parser.addBooleanOption("pairs");

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
        countFlags += (Boolean) parser.getOptionValue(pairedCoverageOpt, false) ? CoverageCounter.PAIRED_COVERAGE : 0;
        String strandopt = (String) parser.getOptionValue(strandOption, "");
        if (strandopt.equals("read")) {
            countFlags += CoverageCounter.STRANDS_BY_READ;
        } else if (strandopt.equals("first")) {
            countFlags += CoverageCounter.STRANDS_BY_FIRST_IN_PAIR;
        } else if (strandopt.equals("second")) {
            log.warn("Warning: 'second' Option undocumented and may be removed in the future. BE WARNED!");
            countFlags += CoverageCounter.STRANDS_BY_SECOND_IN_PAIR;
        }
        return countFlags;
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

        userMessageWriter.println("gct -> igv: " + ifile + " -> " + ofile.getAbsolutePath());

        File tmpDir = null;
        if (tmpDirName != null && tmpDirName.trim().length() > 0) {
            tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
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

        log.info("toTDF.  File = " + ifile);
        log.info("Max zoom = " + maxZoomValue);
        if (probeFile != null && probeFile.trim().length() > 0) {
            log.info("Probe file = " + probeFile);
        }
        String wfString = "Window functions: ";
        for (WindowFunction wf : windowFunctions) {
            wfString += wf.toString() + " ";
        }
        log.info(wfString);

        boolean isGCT = isGCT(typeString);
        Genome genome = loadGenome(genomeId);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }
        File inputFileOrDir = new File(ifile);

        // Estimate the total number of lines to be parsed, for progress updates
        int nLines = estimateLineCount(inputFileOrDir);

        // TODO -- move this block of code out of here, this should be done before calling this method
        // Convert gct files to igv format first
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

            inputFileOrDir = igvFile;
            deleteme = igvFile;
            typeString = ".igv";
        }

        // Convert to tdf
        File outputFile = new File(ofile);
        try {
            Preprocessor p = new Preprocessor(outputFile, genome, windowFunctions, nLines, null);
            if (inputFileOrDir.isDirectory() || inputFileOrDir.getName().endsWith(".list")) {
                p.setSizeEstimate(0);
                List<File> files = getFilesFromDirOrList(inputFileOrDir);
                for (File f : files) {
                    p.preprocess(f, maxZoomValue, typeString);
                }
            } else {
                p.preprocess(inputFileOrDir, maxZoomValue, typeString);
            }
            p.finish();
        } catch (IOException e) {
            log.error(e.getMessage(), e);
            // Delete output file as its probably corrupt
            if (outputFile.exists()) {
                outputFile.delete();
            }
        } finally {
            if (deleteme != null && deleteme.exists()) {
                deleteme.delete();
            }
        }
        userMessageWriter.flush();
    }


    /**
     * Return either (a) the children files in a directory, or (b) files listed in the input, assuming inputFileOrDir
     * is a text file with 1 file listing per line.
     *
     * @param inputFileOrDir
     * @return
     */
    private List<File> getFilesFromDirOrList(File inputFileOrDir) {

        if (inputFileOrDir.isDirectory()) {
            File[] files = inputFileOrDir.listFiles();
            Arrays.sort(files, new Comparator<File>() {
                public int compare(File file, File file1) {
                    return file.getName().compareTo(file1.getName());
                }
            });
            return Arrays.asList(files);
        } else {
            // Must be a "list" file
            BufferedReader br = null;
            try {
                ArrayList<File> files = new ArrayList<File>();
                br = new BufferedReader(new FileReader(inputFileOrDir));
                File parentDirectory = inputFileOrDir.getParentFile();
                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    File f = new File(nextLine);
                    if (f.exists()) {
                        if (f.isDirectory()) {
                            continue; // Skip directories
                        }
                        files.add(f);
                    } else {
                        // Might be relative path
                        f = new File(parentDirectory, nextLine);
                        if (f.exists()) {
                            files.add(f);
                        } else {
                            log.error("File not found: " + nextLine);
                        }
                    }
                }
                return files;
            } catch (IOException e) {
                log.error(e.getMessage(), e);
                throw new RuntimeException("Error parsing input file: " + inputFileOrDir.getAbsolutePath() + " " +
                        e.getMessage());
            } finally {
                if (br != null)
                    try {
                        br.close();
                    } catch (IOException e) {
                        log.error(e.getMessage(), e);
                    }
            }
        }

    }

    private boolean isGCT(String tmp) {
        if (tmp.endsWith(".txt")) tmp = tmp.substring(0, tmp.length() - 4);
        if (tmp.endsWith(".gz")) tmp = tmp.substring(0, tmp.length() - 3);
        return tmp.endsWith(".gct") || tmp.endsWith(".tab") || tmp.equals("mage-tab");
    }

    /**
     * Compute coverage or density of an alignment or feature file.
     *
     * @param ifile           Alignment or feature file
     * @param ofile           Output file
     * @param genomeId        Genome id (e.g. hg18) or full path to a .genome file (e.g. /xchip/igv/scer2.genome)
     * @param maxZoomValue    Maximum zoom level to precompute.  Default value is 7
     * @param windowFunctions
     * @param windowSizeValue
     * @param extFactorValue
     * @param trackLine
     * @param queryString
     * @param minMapQuality
     * @param countFlags
     * @throws IOException
     */
    public void doCount(String ifile, String ofile, String genomeId, int maxZoomValue,
                        Collection<WindowFunction> windowFunctions, int windowSizeValue,
                        int extFactorValue, int preExtFactorValue, int postExtFactorValue,
                        String trackLine, String queryString, int minMapQuality, int countFlags) throws IOException {


        log.info("Computing coverage.  File = " + ifile);
        log.info("Max zoom = " + maxZoomValue);
        log.info("Window size = " + windowSizeValue);
        String wfString = "Window functions: ";
        for (WindowFunction wf : windowFunctions) {
            wfString += wf.toString() + " ";
        }
        log.info(wfString);
        log.info("Ext factor = " + extFactorValue);


        Genome genome = loadGenome(genomeId);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }

        // Multiple files allowed for count command (a tdf and a wig)
        File tdfFile = null;
        File wigFile = null;
        boolean wigStdOut = false;
        String[] files = ofile.split(",");
        
        for(String fileTok: files){
            if (fileTok.endsWith("wig")) {
                wigFile = new File(fileTok);
            } else if (fileTok.endsWith("tdf")) {
                tdfFile = new File(fileTok);
            }else if (fileTok.equals(STDOUT_FILE_STR)){
                wigStdOut = true;
            }
        }

        if (tdfFile != null && !tdfFile.getName().endsWith(".tdf")) {
            tdfFile = new File(tdfFile.getAbsolutePath() + ".tdf");
        }

        try {

            Preprocessor p = new Preprocessor(tdfFile, genome, windowFunctions, -1, null);

            p.setSkipZeroes(true);

            CoverageCounter counter = new CoverageCounter(ifile, p, windowSizeValue, extFactorValue, wigFile,
                    genome, queryString, minMapQuality, countFlags);
            counter.setWriteStdOut(wigStdOut);
            counter.setPreExtFactor(preExtFactorValue);
            counter.setPosExtFactor(postExtFactorValue);

            String prefix = FilenameUtils.getName(ifile);
            String[] tracknames = counter.getTrackNames(prefix + " ");
            p.setTrackParameters(TrackType.COVERAGE, trackLine, tracknames);

            p.setSizeEstimate(((int) (genome.getNominalLength() / windowSizeValue)));

            counter.parse();

            p.finish();

        } catch (Exception e) {
            // Delete the output file(s) as they are probably corrupt
            log.error(e.getMessage(), e);
            if (tdfFile != null && tdfFile.exists()) {
                tdfFile.delete();
            }
            if (tdfFile != null && wigFile.exists()) {
                wigFile.delete();
            }
        }

        userMessageWriter.flush();
    }


    public void doWIBtoWIG(File txtFile, File wibFile, File wigFile, String trackLine) {
        UCSCUtils.convertWIBFile(txtFile, wibFile, wigFile, trackLine);
    }

    public String doIndex(String ifile, String outputDir, int indexType, int binSize) throws IOException {
        String typeString = Preprocessor.getExtension(ifile);
        return doIndex(ifile, typeString, outputDir, indexType, binSize);
    }

    /**
     * Create an index for an alignment or feature file
     * The output index will have the same base name is the input file, although
     * it may be in a different directory. An appropriate index extension (.sai, .idx, etc.) will
     * be appended.
     *
     * @param ifile
     * @param typeString
     * @param outputDir
     * @param indexType
     * @param binSize
     * @throws IOException
     */
    public String doIndex(String ifile, String typeString, String outputDir, int indexType, int binSize) throws IOException {
        File inputFile = new File(ifile);

        if (outputDir == null) {
            outputDir = inputFile.getParent();
        }
        String outputFileName = (new File(outputDir, inputFile.getName())).getAbsolutePath();

        if (typeString.endsWith("gz")) {
            log.error("Cannot index a gzipped file");
            throw new PreprocessingException("Cannot index a gzipped file");
        }

        if (typeString.endsWith("bam")) {
            String msg = "Cannot index a BAM file. Use the samtools package for sorting and indexing BAM files.";
            log.error(msg);
            throw new PreprocessingException(msg);
        }

        String[] fastaTypes = new String[]{"fa", "fna", "fasta"};
        boolean isFasta = false;
        //We have different naming conventions for different index files
        if (typeString.endsWith("sam") && !outputFileName.endsWith(".sai")) {
            outputFileName += ".sai";
        } else if (typeString.endsWith("bam") && !outputFileName.endsWith(".bai")) {
            outputFileName += ".bai";
        } else {

            for (String ft : fastaTypes) {
                if (typeString.endsWith(ft) && !outputFileName.endsWith(".fai")) {
                    outputFileName += ".fai";
                    isFasta = true;
                    break;
                }
            }


            if (!isFasta && !outputFileName.endsWith(".idx")) {
                outputFileName += ".idx";
            }
        }


        File outputFile = new File(outputFileName);

        //Sam/FASTA files are special
        try {
            if (typeString.endsWith("sam")) {
                AlignmentIndexer indexer = AlignmentIndexer.getInstance(inputFile, null, null);
                indexer.createSamIndex(outputFile);
                return outputFileName;
            } else if (isFasta) {
                FastaUtils.createIndexFile(inputFile.getAbsolutePath(), outputFileName);
                return outputFileName;
            }
        } catch (Exception e) {
            // Delete output file as it is probably corrupt
            if (outputFile.exists()) {
                outputFile.delete();
            }
            throw new RuntimeException(e);
        }


        Genome genome = null;  // <= don't do chromosome conversion
        FeatureCodec codec = CodecFactory.getCodec(ifile, genome);
        if (codec != null) {
            try {
                createTribbleIndex(ifile, outputFile, indexType, binSize, codec);
            } catch (TribbleException.MalformedFeatureFile e) {
                StringBuffer buf = new StringBuffer();
                buf.append("<html>Files must be sorted by start position prior to indexing.<br>");
                buf.append(e.getMessage());
                buf.append("<br><br>Note: igvtools can be used to sort the file, select \"File > Run igvtools...\".");
                MessageUtils.showMessage(buf.toString());
            }
        } else {
            throw new DataLoadException("Unknown File Type", ifile);
        }
        userMessageWriter.flush();
        return outputFileName;

    }

    public static void writeTribbleIndex(Index idx, String idxFile) throws IOException{
        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        } catch (Exception e) {
            log.error(e.getMessage(), e);
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
                idxFile = ifile + ".idx";
            }
            writeTribbleIndex(idx, idxFile);
        }
    }


    public void doSort(String ifile, String ofile, String tmpDirName, int maxRecords) {

        userMessageWriter.println("Sorting " + ifile + "  -> " + ofile);
        File inputFile = new File(ifile);
        boolean writeStdOut = ofile.equals(STDOUT_FILE_STR);
        File outputFile = writeStdOut ? null : new File(ofile);
        Sorter sorter = Sorter.getSorter(inputFile, outputFile);
        if (tmpDirName != null && tmpDirName.trim().length() > 0) {
            File tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                log.error("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
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
            if (writeStdOut && outputFile.exists()) {
                outputFile.delete();
            }
        }
        userMessageWriter.flush();
    }


    private void validateArgsLength(String[] nonOptionArgs, int len, String failMessage) throws PreprocessingException {
        if (nonOptionArgs.length < len) {
            throw new PreprocessingException(failMessage + "\n");
        }
    }


    public static Genome loadGenome(String genomeFileOrID) throws IOException {

        String rootDir = FileUtils.getInstallDirectory();

        final GenomeManager genomeManager = GenomeManager.getInstance();
        Genome genome = genomeManager.getCurrentGenome();
        if (genome != null && genome.getId().equals(genomeFileOrID)) {
            return genome;
        }

        File genomeFile = new File(genomeFileOrID);
        if (!genomeFile.exists()) {
            genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID + ".genome");

        }
        if(!genomeFile.exists()) {
            genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID + ".chrom.sizes");
        }
        if (!genomeFile.exists()) {
            genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID);

        }
        if (!genomeFile.exists()) {
            throw new PreprocessingException("Genome definition file not found for: " + genomeFileOrID);
        }


        //TODO Prevents loading genome again if loading from path.
        //May or may not want this, for now we just use it for testing
        if (Globals.isTesting() && genomeFile.getAbsolutePath().endsWith(".genome")) {
            GenomeDescriptor genomeDescriptor = GenomeManager.parseGenomeArchiveFile(genomeFile);
            if (genome != null && genomeDescriptor.getId().equals(genome.getId())) {
                return genome;
            }
        }

        genome = genomeManager.loadGenome(genomeFile.getAbsolutePath(), null);
        if (genome == null) {
            throw new PreprocessingException("Error loading: " + genomeFileOrID);
        }
        return genome;
    }


    /**
     * Test if the file type can be "tiled".
     */
    private static void validateIsTilable(String typeString) {

        if (!(typeString.endsWith("cn") ||
                typeString.endsWith("igv") ||
                typeString.endsWith("wig") ||
                // ifile.toLowerCase().endsWith("cpg.txt") ||
                typeString.endsWith("ewig") ||
                typeString.endsWith("map") ||
                typeString.endsWith("cn") ||
                typeString.endsWith("snp") ||
                typeString.endsWith("xcn") ||
                typeString.endsWith("gct") ||
                typeString.endsWith("tab") ||
                typeString.endsWith("mage-tab") ||
                typeString.endsWith("bedgraph") ||
                typeString.endsWith("ewig.list") ||
                Preprocessor.isAlignmentFile(typeString))) {
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
                if (wf.matches("p\\d{1,2}")) {
                    wf = wf.replaceFirst("p", "percentile");
                }
                try {
                    funcs.add(WindowFunction.valueOf(wf));
                } catch (Exception e) {
                    log.error("Unrecognized window function: " + tokens[i]);
                }
            }
            return funcs;
        }
    }


    /**
     * Estimage the number of lines in the given file, or all files in the given directory, or all files
     * referenced in a ".list" file.
     *
     * @param file a file or directory.
     * @return
     */
    private int estimateLineCount(File file) throws IOException {

        int nLines = 0;
        if (file.isDirectory() || file.getName().endsWith(".list")) {
            List<File> files = getFilesFromDirOrList(file);
            for (File f : files) {
                if (!f.isDirectory()) {
                    nLines += ParsingUtils.estimateLineCount(f.getAbsolutePath());
                }
            }
        } else if (file.getName().endsWith(".map")) {
            BufferedReader reader = null;
            try {
                reader = new BufferedReader(new FileReader(file));
                String nextLine;

                while ((nextLine = reader.readLine()) != null) {
                    String[] tokens = Globals.whitespacePattern.split(nextLine);
                    File f = new File(file.getParent(), tokens[1]);
                    nLines += ParsingUtils.estimateLineCount(f.getAbsolutePath());
                }
                return nLines;

            } catch (Exception e) {
                userMessageWriter.println("Error estimating line count: " + e.getMessage());
                nLines = 0;
            } finally {
                if (reader != null) reader.close();
            }
        } else {
            nLines = ParsingUtils.estimateLineCount(file.getAbsolutePath());
        }
        return nLines;

    }


    private static void launchGUI() {
        IgvToolsGui.main(null);
    }

}

