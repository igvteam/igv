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
package org.broad.igv.tools;


import jargs.gnu.CmdLineParser;
import org.broad.igv.Globals;
import org.broad.igv.data.seg.SegmentedDataWriter;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.reader.AlignmentIndexer;
import org.broad.igv.tools.converters.ExpressionFormatter;
import org.broad.igv.tools.converters.GCTtoIGVConverter;
import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.track.tribble.CodecFactory;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.vcf.VCFtoBed;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Command line parser for "igvtools".
 *
 * @author jrobinso
 */
public class IgvTools {

    static String version = "@VERSION";

    static String[] commandDocs = new String[]{
            "version print the version number",
            "sort    sort an alignment file by start position",
            "index   index an alignment file",
            "tile    convert an input file (cn, gct, wig) to tiled data format (tdf)",
            "count   compute coverage density for an alignment file"
    };
    static final int MAX_RECORDS_IN_RAM = 500000;
    public static final int MAX_ZOOM = 7;
    public static final int WINDOW_SIZE = 25;
    public static final int EXT_FACTOR = 0;
    public static final Object PROBE_FILE = null;
    public static final int LINEAR_BIN_SIZE = 16000;
    public static final int INTERVAL_SIZE = 1000;
    public static final int LINEAR_INDEX = 1;
    public static final int INTERVAL_INDEX = 2;

    /**
     * The general usage string
     */
    static String usageString() {
        StringBuffer buf = new StringBuffer();

        buf.append("\nProgram: igvtools\n\n");
        buf.append("Usage: igvtools [command] [options] [arguments]\n\n");
        buf.append("Command:");
        for (String c : commandDocs) {
            buf.append(" " + c + "\n\t");

        }
        return buf.toString();
    }


    /**
     * Main method.  Used to run igvtools from the command line
     *
     * @param argv
     * @throws IOException
     * @throws PreprocessingException
     */
    public static void main(String[] argv) throws IOException, PreprocessingException {

        Globals.setHeadless(true);

        (new IgvTools()).run(argv);
    }

    private void parseArguments(String[] argv) {

    }

    private void run(String[] argv) {
        CmdLineParser parser = new CmdLineParser();
        CmdLineParser.Option helpOption = parser.addBooleanOption('h', "help");
        CmdLineParser.Option guiOption = parser.addBooleanOption('g', "gui");

        // general options
        CmdLineParser.Option windowFunctions = parser.addStringOption('f', "windowFunctions");
        CmdLineParser.Option tmpDirOption = parser.addStringOption('t', "tmpDir");
        CmdLineParser.Option maxZoomOption = parser.addIntegerOption('z', "maxZoom");
        CmdLineParser.Option compressOption = parser.addBooleanOption('c', "compressed");

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

        // options for gentest
        //CmdLineParser.Option nRows = parser.addIntegerOption('r', "numRows");
        CmdLineParser.Option strandOption = parser.addIntegerOption('s', "strand");
        //CmdLineParser.Option unOrdered = parser.addBooleanOption('u', "unOrdered");

        // options for coverage
        CmdLineParser.Option coverageOptions = parser.addStringOption('x', "coverageOptions");


        // Parse optional arguments (switches, etc)
        try {
            parser.parse(argv);
        } catch (CmdLineParser.OptionException e) {
            System.err.println(e.getMessage());
            return;
        }

        boolean help = ((Boolean) parser.getOptionValue(helpOption, false));
        if (argv.length < 1 || help) {
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
            validateArgsLength(nonOptionArgs, 1);

            // The command
            String command = nonOptionArgs[EXT_FACTOR].toLowerCase();

            // Do "version" now, its the only command with no arguments
            if (command.equals("version")) {
                System.out.println("Version " + version);
                return;
            }

            // All remaining commands require an input file, and most need the file extension.  Do that here.
            validateArgsLength(nonOptionArgs, 2);
            String ifile = nonOptionArgs[1];
            if (!(new File(ifile)).exists()) {
                throw new PreprocessingException("File not found: " + ifile);
            }

            int maxRecords = (Integer) parser.getOptionValue(maxRecordsOption, MAX_RECORDS_IN_RAM);

            if (command.equals("count") || command.equals("pre") || command.equals("tile")) {
                validateArgsLength(nonOptionArgs, 4);
                int windowSizeValue = (Integer) parser.getOptionValue(windowSizeOption, WINDOW_SIZE);
                int maxZoomValue = (Integer) parser.getOptionValue(maxZoomOption, MAX_ZOOM);
                String ofile = nonOptionArgs[2];
                String genomeId = nonOptionArgs[3];
                boolean isGCT = Preprocessor.getExtension(ifile).endsWith("gct");
                String wfsString = (String) parser.getOptionValue(windowFunctions);
                Collection<WindowFunction> wfList = parseWFS(wfsString, isGCT);
                String coverageOpt = (String) parser.getOptionValue(coverageOptions);

                if (command.equals("count")) {
                    int extFactorValue = (Integer) parser.getOptionValue(extFactorOption, EXT_FACTOR);
                    int strandOptionValue = (Integer) parser.getOptionValue(strandOption, -1);
                    doCount(ifile, ofile, genomeId, maxZoomValue, wfList, windowSizeValue, extFactorValue,
                            strandOptionValue, coverageOpt);
                } else {
                    String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                    doTile(ifile, ofile, probeFile, genomeId, maxZoomValue, wfList, tmpDirName, maxRecords);
                }
            } else if (command.toLowerCase().equals("sort")) {
                validateArgsLength(nonOptionArgs, 3);

                String ofile = nonOptionArgs[2];
                doSort(ifile, ofile, tmpDirName, maxRecords);
            } else if (command.equals("index")) {
                int indexType = (Integer) parser.getOptionValue(indexTypeOption, LINEAR_INDEX);
                int defaultBinSize = indexType == LINEAR_INDEX ? LINEAR_BIN_SIZE : INTERVAL_SIZE;
                int binSize = (Integer) parser.getOptionValue(binSizeOption, defaultBinSize);
                String outputDir = (String) parser.getOptionValue(outputDirOption, null);
                doIndex(ifile, outputDir, indexType, binSize);
            } else if (command.equals("seg")) {
                doPreprocessSeg(argv);
            } else if (command.equals("wibtowig")) {
                validateArgsLength(nonOptionArgs, 4);
                File txtFile = new File(nonOptionArgs[1]);
                File wibFile = new File(nonOptionArgs[2]);
                File wigFile = new File(nonOptionArgs[3]);
                String trackLine = nonOptionArgs.length > 4 ? nonOptionArgs[4] : null;
                doWIBtoWIG(txtFile, wibFile, wigFile, trackLine);
            } else if (command.equals("splitgff")) {
                validateArgsLength(nonOptionArgs, 3);
                String outputDirectory = nonOptionArgs[2];
                GFFParser.splitFileByType(ifile, outputDirectory);
                //} else if (command.equals("gentest")) {
                //  validateArgsLength(nonOptionArgs, 2);
                //  int nRowsValue = (Integer) parser.getOptionValue(nRows, 10);
                //  boolean sorted = !((Boolean) parser.getOptionValue(unOrdered, false));
                //  int strandOptionValue = (Integer) parser.getOptionValue(strandOption, -1);
                //  String ofile = ifile;
                //  TestFileGenerator.generateTestFile(ofile, sorted, nRowsValue, strandOptionValue);
            } else if (command.toLowerCase().equals("gcttoigv")) {
                validateArgsLength(nonOptionArgs, 4);
                String ofile = nonOptionArgs[2];
                // Output files must have .igv extension
                if (!ofile.endsWith(".igv")) {
                    ofile = ofile + ".igv";
                }
                String genomeId = nonOptionArgs[3];
                String probeFile = (String) parser.getOptionValue(probeFileOption, PROBE_FILE);
                doGCTtoIGV(ifile, new File(ofile), probeFile, genomeId, maxRecords, tmpDirName);
            } else if (command.equals("formatexp")) {
                validateArgsLength(nonOptionArgs, 2);
                File inputFile = new File(nonOptionArgs[1]);
                File outputFile = new File(nonOptionArgs[2]);
                (new ExpressionFormatter()).convert(inputFile, outputFile);
            } else if (command.equals("wigtobed")) {
                validateArgsLength(nonOptionArgs, 1);
                String inputFile = nonOptionArgs[1];
                float hetThreshold = 0.17f;
                if(nonOptionArgs.length > 2) {
                    hetThreshold = Float.parseFloat(nonOptionArgs[2]);
                }
                float homThreshold = 0.55f;
                if(nonOptionArgs.length > 3) {
                    homThreshold = Float.parseFloat(nonOptionArgs[3]);
                }
                WigToBed.run(inputFile, hetThreshold, homThreshold);
            } else if (command.equals("vcftobed")) {
                validateArgsLength(nonOptionArgs, 2);
                String inputFile = nonOptionArgs[1];
                String outputFile = nonOptionArgs[2];
                VCFtoBed.convert(inputFile, outputFile);
            } else if(command.equals("lanecounter")) {
                validateArgsLength(nonOptionArgs, 3);
                Genome genome = loadGenome(nonOptionArgs[1]);
                String bamFileList = nonOptionArgs[2];
                String queryInterval = nonOptionArgs[3];
                LaneCounter.run(genome, bamFileList, queryInterval);
            }
            else {
                throw new PreprocessingException("Unknown command: " + argv[EXT_FACTOR]);
            }
        } catch (PreprocessingException e) {
            System.err.println(e.getMessage());
        } catch (IOException e) {
            throw new PreprocessingException("Unexpected IO error: ", e);
        }
    }

    private void doGCTtoIGV(String ifile, File ofile, String probefile, String genomeId, int maxRecords, String tmpDirName) throws IOException {
        Genome genome = loadGenome(genomeId);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }

        File tmpDir = null;
        if (tmpDirName != null && tmpDirName.trim().length() > 0) {
            tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                System.err.println("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
                throw new PreprocessingException("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
            }
        }

        GCTtoIGVConverter.convert(new File(ifile), ofile, probefile, genomeId, maxRecords, tmpDir);

    }

    public void doTile(String ifile, String ofile, String probeFile, String genomeId, int maxZoomValue,
                       Collection<WindowFunction> windowFunctions) throws IOException, PreprocessingException {
        doTile(ifile, ofile, probeFile, genomeId, maxZoomValue, windowFunctions, null, MAX_RECORDS_IN_RAM);
    }


    public void doTile(String ifile, String ofile, String probeFile, String genomeId, int maxZoomValue,
                       Collection<WindowFunction> windowFunctions, String tmpDirName, int maxRecords) throws IOException, PreprocessingException {
        validateIsTilable(ifile);

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

        Genome genome = loadGenome(genomeId);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }
        File tmp = new File(ifile);


        int nLines = tmp.isDirectory() ? ParsingUtils.estimateLineCount(tmp) : ParsingUtils.estimateLineCount(ifile);

        // Convert  gct files to igv format first
        if (isGCT(ifile)) {
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
            doGCTtoIGV(ifile, igvFile, probeFile, genomeId, maxRecords, tmpDirName);

            tmp = igvFile;

        }


        Preprocessor p = new Preprocessor(new File(ofile), genome, windowFunctions, nLines, null);
        if (tmp.isDirectory()) {
            for (File f : tmp.listFiles()) {
                p.preprocess(f, probeFile, maxZoomValue);
            }
        } else {
            p.preprocess(tmp, probeFile, maxZoomValue);
        }
        p.finish();

        System.out.flush();

    }

    private boolean isGCT(String iFile) {
        String tmp = iFile;
        if (tmp.endsWith(".txt")) tmp = tmp.substring(0, tmp.length() - 4);
        if (tmp.endsWith(".gz")) tmp = tmp.substring(0, tmp.length() - 3);
        return tmp.endsWith(".gct") || tmp.endsWith(".tab");
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
     * @param strandOption
     * @throws IOException
     */
    public void doCount(String ifile, String ofile, String genomeId, int maxZoomValue,
                        Collection<WindowFunction> windowFunctions,
                        int windowSizeValue, int extFactorValue, int strandOption,
                        String coverageOpt) throws IOException {


        System.out.println("Computing coverage.  File = " + ifile);
        System.out.println("Max zoom = " + maxZoomValue);
        System.out.println("Window size = " + windowSizeValue);
        System.out.print("Window functions: ");
        for (WindowFunction wf : windowFunctions) {
            System.out.print(wf.toString() + " ");
        }
        System.out.println();
        System.out.println("Ext factor = " + extFactorValue);


        Genome genome = loadGenome(genomeId);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: " + genomeId);
        }
        int nLines = ParsingUtils.estimateLineCount(ifile);

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

        Preprocessor p = new Preprocessor(tdfFile, genome, windowFunctions, nLines, null);
        p.count(ifile, windowSizeValue, extFactorValue, maxZoomValue, wigFile, strandOption, coverageOpt);
        p.finish();

        System.out.flush();
    }


    public void doWIBtoWIG(File txtFile, File wibFile, File wigFile, String trackLine) {
        UCSCUtils.convertWIBFile(txtFile, wibFile, wigFile, trackLine);
    }


    /**
     * @param argv
     * @deprecated Convert a .seg file to a .seg.zip file
     */
    public void doPreprocessSeg(String[] argv) {
        int nArgs = argv.length - 1;
        String[] args = new String[nArgs];
        System.arraycopy(argv, 1, args, EXT_FACTOR, nArgs);
        SegmentedDataWriter.main(args);
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
        if (ifile.endsWith(".bed") || ifile.endsWith(".vcf") || ifile.endsWith(".vcf4") || ifile.endsWith(".gff") ||
                ifile.endsWith(".gff3") || ifile.endsWith(".psl") || ifile.endsWith(".pslx")) {

            if (!outputFileName.endsWith(".idx")) {
                outputFileName = outputFileName + ".idx";
            }
            try {
                createTribbleIndex(ifile, new File(outputFileName), indexType, binSize);
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
            indexer.createSamIndex(new File(outputFileName));
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
    private void createTribbleIndex(String ifile, File outputFile, int indexType, int binSize) throws IOException {
        FeatureCodec codec = CodecFactory.getCodec(ifile);
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
            }
            finally {
                if (stream != null) {
                    stream.close();
                }
            }
        }
    }


    public void doSort(String ifile, String ofile, String tmpDirName, int maxRecords) {

        System.out.println("Sorting " + ifile);
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
        sorter.run();
        System.out.println("Done");
        System.out.flush();
    }


    private void validateArgsLength(String[] nonOptionArgs, int len) throws PreprocessingException {
        if (nonOptionArgs.length < len) {
            throw new PreprocessingException(usageString());
        }
    }


    public static Genome loadGenome(String genomeFileOrID) {


        String rootDir = FileUtils.getInstallDirectory();

        Genome genome = GenomeManager.getInstance().getGenome(genomeFileOrID);
        if (genome != null) {
            return genome;
        }


        File genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID + ".genome");
        if (!genomeFile.exists()) {
            genomeFile = new File(rootDir, "genomes" + File.separator + genomeFileOrID);
        }
        if (!genomeFile.exists()) {
            genomeFile = new File(genomeFileOrID);
        }
        if (!genomeFile.exists()) {
            throw new PreprocessingException("Genome definition file not found for: " + genomeFileOrID);
        }

        GenomeManager.GenomeListItem item = GenomeManager.getInstance().loadGenomeFromLocalFile(genomeFile);
        String genomeId = item.getId();
        if (genomeId == null) {
            throw new PreprocessingException("Error loading: " + genomeFileOrID);
        }
        GenomeManager.getInstance().setGenomeId(genomeId);
        genome = GenomeManager.getInstance().getGenome(genomeId);
        return genome;
    }


    /**
     * Test if the file type can be "tiled".
     */
    private static void validateIsTilable(String ifile) {
        File tmp = new File(ifile);
        if (tmp.isDirectory()) {
            File[] files = tmp.listFiles();
            if (files.length == EXT_FACTOR) {
                throw new PreprocessingException("Tile command not supported for empty directory: " + ifile);
            }
            ifile = files[EXT_FACTOR].getName();
        }
        String ext = Preprocessor.getExtension(ifile);
        if (!(ext.equals(".cn") ||
                ext.equals(".igv") ||
                ext.equals(".wig") ||
                ext.equals(".ewig") ||
                ext.equals(".cn") ||
                ext.equals(".snp") ||
                ext.equals(".xcn") ||
                ext.equals(".gct") ||
                ext.equals(".bedgraph") ||
                Preprocessor.isAlignmentFile(ext))) {
            throw new PreprocessingException("Tile command not supported for files of type: " + ext);
        }
    }


    /**
     * Parse the window functions line.
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

