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

package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.Dataset;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.data.WiggleParser;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.sam.reader.SamUtils;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.igv.tools.sort.SorterTest;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.junit.*;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.*;
import java.util.*;

import static junit.framework.Assert.*;

public class IGVToolsTest extends AbstractHeadlessTest {

    IgvTools igvTools;

    private static final String hg18id = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";
    private static final int MAX_LINES_CHECK = 200;

    @Rule
    public TestRule testTimeout = new Timeout((int) 1e3*60*60);

    @Before
    public void setUp() throws Exception {
        igvTools = new IgvTools();

    }

    @After
    public void tearDown() throws Exception {
        igvTools = null;
    }


    @Test
    public void testIndexSam() throws Exception {
        String samFile = TestUtils.DATA_DIR + "sam/NA12878.muc1.test2.sam";
        String indPath = samFile + ".sai";
        File indFile = new File(indPath);
        indFile.delete();
        indFile.deleteOnExit();


        igvTools.doIndex(samFile, indPath, IgvTools.LINEAR_INDEX, IgvTools.LINEAR_BIN_SIZE);

        //Check that only the index file we intended exists
        assertTrue(indFile.exists());
        assertFalse((new File(samFile + ".idx").exists()));
        assertFalse((new File(indPath + ".idx").exists()));

        FeatureIndex idx = SamUtils.getIndexFor(samFile);
        assertTrue(idx.containsChromosome("chr1"));
        assertEquals(1, idx.getIndexedChromosomes().size());
    }

    @Test
    public void testLinearIndex() throws IOException {

        String bedFile = TestUtils.DATA_DIR + "bed/test.bed";

        File idxFile = new File(bedFile + ".idx");
        if (idxFile.exists()) {
            idxFile.delete();
        }

        igvTools.doIndex(bedFile, 1, 16000);

        assertTrue(idxFile.exists());

        Index idx = IndexFactory.loadIndex(idxFile.getAbsolutePath());

        List<Block> blocks = idx.getBlocks("chr1", 100, 200);
        Block block = blocks.get(0);
        assertEquals("Unexpected start position ", 0, block.getStartPosition());

    }

    @Test
    public void testIntervalIndex33() throws Exception {
        String testFile = TestUtils.LARGE_DATA_DIR + "CEU.SRP000032.2010_03_v3.3.genotypes.head.vcf";
        FeatureCodec codec = new VCF3Codec();
        tstIntervalIndex(testFile, codec);
    }

    @Test
    public void testIntervalIndex40() throws Exception {
        String testFile = TestUtils.LARGE_DATA_DIR + "CEU.SRP000032.2010_03_v4.0.genotypes.head.vcf";
        FeatureCodec codec = new VCFCodec();
        tstIntervalIndex(testFile, codec);
    }

    private void tstIntervalIndex(String testFile, FeatureCodec codec) throws IOException {

        // Create an interval tree index with 5 features per interval
        File indexFile = new File(testFile + ".idx");
        if (indexFile.exists()) {
            indexFile.delete();
        }
        igvTools.doIndex(testFile, 2, 5);
        indexFile.deleteOnExit();

        // Now use the index
        String chr = "1";
        int start = 1718546;
        int end = 1748915;
        int[] expectedStarts = {1718547, 1718829, 1723079, 1724830, 1731376, 1733967, 1735586, 1736016, 1738594,
                1739272, 1741124, 1742815, 1743224, 1748886, 1748914};

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(testFile, codec);
        Iterator<VariantContext> iter = bfr.query(chr, start, end);
        int count = 0;
        while (iter.hasNext()) {
            VariantContext feat = iter.next();
            int expStart = expectedStarts[count];
            assertEquals(expStart, feat.getStart());
            count++;
        }
        Assert.assertEquals(15, count);
    }


    @Test
    public void testVersion() throws IOException {
        String[] args = {"version"};
        //IgvTools.main(args);
        igvTools.run(args);
    }

    @Test
    public void testTileWigFile() throws IOException {

        String inputFile = TestUtils.LARGE_DATA_DIR + "phastCons_chr1.wig";
        testTile(inputFile, 0, 0);
    }

    @Test
    public void testTileCNFile() throws IOException {

        String inputFile = TestUtils.DATA_DIR + "cn/HindForGISTIC.hg16.cn";
        testTile(inputFile, 5000000, 6000000);
    }


    @Test
    public void testTileGCT() throws IOException {
        String inputFile = TestUtils.DATA_DIR + "gct/OV.transcriptome__agilentg4502.data.txt";
        String outFilePath = TestUtils.DATA_DIR + "out/testTileGCT.wig";
        String genome = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";
        String[] args = {"tile", "-z", "1", "--fileType", "mage-tab", inputFile, outFilePath, hg18id};
        igvTools.run(args);

        inputFile = TestUtils.DATA_DIR + "gct/GBM.methylation__sampled.data.txt";
        args = new String[]{"tile", "-z", "1", "--fileType", "mage-tab", inputFile, outFilePath, hg18id};
        igvTools.run(args);

    }


    private void testTile(String inputFile, int start, int end) throws IOException {
        String file1 = TestUtils.DATA_DIR + "out/file1.tdf";
        String file2 = TestUtils.DATA_DIR + "out/file2.tdf";

        //todo Compare 2 outputs more meaningfully
        String[] args = {"toTDF", "-z", "1", "--windowFunctions", "min", inputFile, file1, hg18id};
        igvTools.run(args);

        args = new String[]{"toTDF", "-z", "2", "--windowFunctions", "max", inputFile, file2, hg18id};
        (new IgvTools()).run(args);


        String dsName = "/chr1/raw";

        TDFDataset ds1 = TDFReader.getReader(file1).getDataset(dsName);
        TDFDataset ds2 = TDFReader.getReader(file2).getDataset(dsName);

        TDFTile t1 = ds1.getTiles(start, end).get(0);
        TDFTile t2 = ds2.getTiles(start, end).get(0);

        int nPts = t1.getSize();
        assertEquals(nPts, t2.getSize());

        for (int i = 0; i < nPts; i++) {
            assertTrue(t1.getStartPosition(i) < t1.getEndPosition(i));
            assertEquals(t1.getStartPosition(i), t2.getStartPosition(i));
            assertTrue(t1.getValue(0, i) <= t2.getValue(0, i));
            if (i < nPts - 1) {
                assertTrue(t1.getStartPosition(i) < t1.getStartPosition(i + 1));
            }
        }

        (new File(file1)).delete();
        (new File(file2)).delete();
    }

    @Test
    public void testCountTDF() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        tstCount(inputFile, "testtdf", "tdf", null, -1, -1);
    }

    @Test
    public void testCountWIG() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        tstCount(inputFile, "testwig", "wig", null, -1, -1);
        tstCount(inputFile, "testwig", "wig", "chr2", 178709699, 179008373);
    }

    @Test
    public void testCountSAM() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "sam/test_2.sam";
        tstCount(inputFile, "testwig", "wig", null, -1, -1);
    }

    @Test
    public void testCountBAM() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bam/NA12878.SLX.sample.bam";
        tstCount(inputFile, "testwig", "wig", null, -1, -1);
    }

    public void tstCount(String inputFile, String outputBase, String outputExt,
                         String chr, int start, int end) throws Exception {
        String outputFile = TestUtils.DATA_DIR + "out/" + outputBase + "_";

        boolean query = chr != null && start >= 0 && end >= start + 1;


        String[] opts = new String[]{"", "--bases", "--strands=read", "--strands=first", "--strands=read --bases"};
        Map<String, float[]> rowTotals = new HashMap<String, float[]>(opts.length);
        String refOpt = null;

        for (int ind = 0; ind < opts.length; ind++) {
            String opt = opts[ind];
            int winsize = 5;
            if (query) {
                opt += " --windowSize " + winsize + " --query " + chr + ":" + start + "-" + end;
            }

            if (ind == 0) {
                refOpt = opt;
            }

            String fullout = outputFile + ind + "." + outputExt;
            String input = "count " + opt + " " + inputFile + " " + fullout + " " + hg18id;
            String[] args = input.split("\\s+");
            igvTools.run(args);

            if (outputExt.equals("tdf")) {
                TDFReader reader = TDFReader.getReader(fullout);
                assertTrue(reader.getDatasetNames().size() > 0);
                if (query) {
                    for (String name : reader.getDatasetNames()) {
                        TDFDataset ds = reader.getDataset(name);
                        List<TDFTile> tiles = ds.getTiles();
                        for (TDFTile tile : tiles) {
                            assertTrue(tile.getTileStart() >= start);
                            assertTrue(tile.getTileStart() < end);
                        }
                    }
                }
            } else {
                File outFile = new File(fullout);
                assertTrue(outFile.exists());
                assertTrue(outFile.canRead());

                ResourceLocator locator = new ResourceLocator(fullout);
                WiggleDataset ds = (new WiggleParser(locator, IgvTools.loadGenome(hg18id, true))).parse();

                //We miss a few alignments with this option sometimes,
                //so it doesn't total up the same
                if (!opt.contains("strands=first")) {
                    float[] linetotals = getLineTotals(fullout);
                    rowTotals.put(opt, linetotals);
                }

                if (query) {
                    assertEquals(1, ds.getChromosomes().length);
                    assertEquals(chr, ds.getChromosomes()[0]);

                    int[] starts = ds.getStartLocations(chr);
                    for (Integer act_start : starts) {
                        assertTrue(act_start + " is outside range", act_start >= start - winsize && act_start < end + winsize);
                    }
                }
            }
        }

        //Compare row totals
        //They should all add up the same
        float[] refTotals = rowTotals.get(refOpt);
        for (String opt : rowTotals.keySet()) {
            if (opt.equals(refOpt)) {
                continue;
            }
            float[] testTotals = rowTotals.get(opt);
            assertEquals(refTotals.length, testTotals.length);
            for (int rw = 0; rw < refTotals.length; rw++) {
                float diff = refTotals[rw] - testTotals[rw];
                String msg = "Difference between " + refTotals[rw] + " and " + testTotals[rw] + " too high. ";
                msg += "Row " + rw + ". Test option " + opt;
                //This can get pretty high, we only use 2 digits of precision
                //Also test this in CoverageCounterTest
                assertTrue(msg, Math.abs(diff) <= 0.1);
            }
        }
    }

    /**
     * Calculates the sum of each row, excluding the first column.
     * Skips non-numeric rows
     *
     * @param filename
     * @return
     */
    private float[] getLineTotals(String filename) throws Exception {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";

        float tmpsum;
        List<Float> sums = new ArrayList<Float>();
        while ((line = reader.readLine()) != null && sums.size() < MAX_LINES_CHECK) {
            try {
                String[] tokens = line.split("\\t");
                tmpsum = 0;
                for (int ii = 1; ii < tokens.length; ii++) {
                    tmpsum += Float.parseFloat(tokens[ii]);
                }
                sums.add(tmpsum);
            } catch (NumberFormatException e) {
                continue;
            }

        }

        reader.close();

        float[] toret = new float[sums.size()];
        for (int ii = 0; ii < sums.size(); ii++) {
            toret[ii] = sums.get(ii);
        }

        return toret;

    }

    public static String[] generateRepLargebamsList(String listPath, String bamFiName, int reps) throws IOException {
        return generateRepLargebamsList(listPath, bamFiName, reps, false);
    }

    /*
    Generate a bam.list file, with rep entries, all having the same content bamPath.
     If makeAbsolute is true and bamPath not absolute, the listPath parent directory is prepended.
     The file is saved to listPath
     */
    public static String[] generateRepLargebamsList(String listPath, String bamPath, int reps, boolean makeAbsolute) throws IOException {

        File listFile = new File(listPath);
        listFile.delete();
        listFile.deleteOnExit();
        File f = new File(bamPath);
        String eachPath = null;
        if (makeAbsolute && !f.isAbsolute()) {
            eachPath = FileUtils.getAbsolutePath(bamPath, listPath);
        } else {
            eachPath = f.getPath();
        }
        //We generate the file on each test, because largedata dir can change
        List<String> largebams = new ArrayList<String>(reps);
        for (int ii = 0; ii < reps; ii++) {
            largebams.add(eachPath);
        }
        FileWriter writer = new FileWriter(listFile);
        for (String s : largebams) {
            writer.write(s + "\n");
        }
        writer.close();

        return largebams.toArray(new String[0]);
    }

    /**
     * Test counting from a BAM.list file. Note that if the paths
     * in the bam.list are relative, they are interpreted as relative to
     * the location of the bam.list file, rather than the working directory.
     * <p/>
     * If paths are entered on the command line (option to be deprecated shortly)
     * via comma separated list, they are interpreted relative to the current working directory.
     *
     * @throws Exception
     */
    @Test
    public void testCountBAMList() throws Exception {
        String listPath = TestUtils.DATA_DIR + "bam/test.unindexed.bam.list";
        tstCountBamList(listPath);
    }

    /**
     * Test iterating through a merged bam file(actually a list with the same file duplicated).  Each record
     * should appear twice (since the list contains the same file twice), and in coordinate sort order.
     *
     * @throws Exception
     */
    @Test
    public void testMergedBam() throws Exception {
        String listPath = TestUtils.DATA_DIR + "bam/test.unindexed.bam.list";
        AlignmentReader reader = AlignmentReaderFactory.getReader(new ResourceLocator(listPath), false);

        Set<String> visitedChromosomes = new HashSet();
        String lastChr = null;
        int lastStart = -1;

        Iterator<Alignment> iter = reader.iterator();
        while (iter.hasNext()) {
            Alignment a1 = iter.next();
            Alignment a2 = iter.next();
            assertEquals(a1.getReadName(), a2.getReadName());
            assertEquals(a1.getChr(), a2.getChr());
            assertEquals(a1.getStart(), a2.getStart());
            assertEquals(a2.getEnd(), a2.getEnd());

            String chr = a1.getChr();
            int start = a1.getAlignmentStart();
            if (lastChr != null && chr.equals(lastChr)) {
                assertTrue(a1.getReadName(), start >= lastStart);
            } else {
                assertFalse(visitedChromosomes.contains(chr));
                lastChr = chr;
                visitedChromosomes.add(chr);
            }
            lastStart = start;

        }
    }

    private void tstCountBamList(String listArg) throws Exception {
        String outputFile = TestUtils.DATA_DIR + "out/file_";

        String[] opts = new String[]{"--strands=read", "--strands=first", ""};

        for (int ind = 0; ind < opts.length; ind++) {
            String opt = opts[ind];
            String fullout = outputFile + ind + ".tdf";
            String input = "count " + opt + " " + listArg + " " + fullout + " " + hg18id;
            String[] args = input.split("\\s+");
            igvTools.run(args);

            TDFReader reader = TDFReader.getReader(fullout);
            assertTrue(reader.getDatasetNames().size() > 0);
        }
    }

    @Test
    public void testSort() throws Exception {
        String inputFiname = "Unigene.unsorted.bed";
        String inputFile = TestUtils.DATA_DIR + "bed/" + inputFiname;
        String outputFile = TestUtils.DATA_DIR + "out/" + inputFiname + ".sorted";
        File oFile = new File(outputFile);
        oFile.deleteOnExit();

        String input = "sort --tmpDir=./ --maxRecords=50 " + inputFile + " " + outputFile;
        igvTools.run(input.split("\\s+"));
        int numlines = SorterTest.checkBedSorted(oFile);
        assertEquals(71, numlines);
    }

    /**
     * This test could stand to be improved, but it's difficult to test math.
     * So we just check that file is about the right size (and well formed).
     *
     * @throws Exception
     */
    @Test
    public void testFormatexp() throws Exception {
        String inputFiname = "igv_test2";
        String ext = ".gct";
        String inputFile = TestUtils.DATA_DIR + "gct/" + inputFiname + ext;
        String outputFile = TestUtils.DATA_DIR + "out/" + inputFiname + "_formatted" + ext;
        File oFile = new File(outputFile);
        oFile.deleteOnExit();

        String input = "formatexp " + inputFile + " " + outputFile;
        igvTools.run(input.split("\\s+"));
        Genome genome = TestUtils.loadGenome();

        ExpressionFileParser parser = new ExpressionFileParser(new ResourceLocator(outputFile), null, genome);
        Dataset ds = parser.createDataset();
        assertEquals(10, ds.getChromosomes().length);
    }

    @Test
    public void testCountDups() throws Exception {
        String inputFiname = "test_5duplicates";
        String ext = ".sam";
        String inputFile = TestUtils.DATA_DIR + "sam/" + inputFiname + ext;
        String outputFileND = TestUtils.DATA_DIR + "out/" + inputFiname + "_nodups" + ".tdf";
        String outputFileWithDup = TestUtils.DATA_DIR + "out/" + inputFiname + "_withdups" + ".tdf";

        String queryChr = "1";

        int pos = 9718611;
        String queryStr = queryChr + ":" + (pos - 100) + "-" + (pos + 100) + " ";
        String cmd_nodups = "count --windowSize 1 -z 7 --query " + queryStr + inputFile + " " + outputFileND + " " + hg18id;
        igvTools.run(cmd_nodups.split("\\s+"));

        String cmd_withdups = "count --includeDuplicates -z 7 --windowSize 1 --query " + queryStr + inputFile + " " + outputFileWithDup + " " + hg18id;
        igvTools.run(cmd_withdups.split("\\s+"));

        assertTrue((new File(outputFileND).exists()));
        assertTrue((new File(outputFileWithDup).exists()));

        Genome genome = TestUtils.loadGenome();

        //Have to read back in using aliased chromosome names
        String readChr = genome.getChromosomeAlias(queryChr);

        int noDupCount = (int) getCount(outputFileND, readChr, 23, pos, genome);
        int dupCount = (int) getCount(outputFileWithDup, readChr, 23, pos, genome);

        assertEquals(noDupCount + 4, dupCount);

        //No dups at this location
        pos += 80;
        noDupCount = (int) getCount(outputFileND, readChr, 23, pos, genome);
        dupCount = (int) getCount(outputFileWithDup, readChr, 23, pos, genome);
        assertEquals(noDupCount, dupCount);

    }

    private float getCount(String filename, String chr, int zoom, int pos, Genome genome) {
        TDFReader reader = TDFReader.getReader(filename);
        TDFDataset ds = reader.getDataset(chr, zoom, WindowFunction.mean);
        TDFDataSource dataSource = new TDFDataSource(reader, 0, "test", genome);
        List<LocusScore> scores = dataSource.getSummaryScoresForRange(chr, pos - 1, pos + 1, zoom);
        return scores.get(0).getScore();
    }

    @Test
    public void testCountAliased() throws Exception {
        tstCountAliased(TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam",
                TestUtils.DATA_DIR + "sam/NA12878.muc1.test_modchr.sam");
    }


    /**
     * Test count when using a custom alias file.
     * Should supply a file using normal chromosome names, and aliases
     * which are defined in the genome file.
     */
    public void tstCountAliased(String normfile, String aliasedfile) throws Exception {
        String genfile = TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome";
        String outfile = TestUtils.DATA_DIR + "out/tmpcount1.wig";
        File outFi = new File(outfile);
        outFi.delete();
        outFi.deleteOnExit();

        //Count aliased file
        String command = "count " + normfile + " " + outfile + " " + genfile;
        igvTools.run(command.split("\\s+"));

        //Count non-aliased file
        String outfile2 = TestUtils.DATA_DIR + "out/tmpcount2.wig";
        File outFi2 = new File(outfile2);
        outFi2.delete();
        //Count aliased file
        command = "count " + aliasedfile + " " + outfile2 + " " + genfile;

        igvTools = new IgvTools();
        igvTools.run(command.split("\\s+"));

        BufferedReader reader1 = new BufferedReader(new FileReader(outfile));
        BufferedReader reader2 = new BufferedReader(new FileReader(outfile2));
        String line1, line2;
        line1 = reader1.readLine();
        line2 = reader2.readLine();
        while ((line1 = reader1.readLine()) != null) {
            line2 = reader2.readLine();
            //Header lines don't need to match
            if (line1.startsWith("#") || line1.startsWith("variableStep")
                    || line1.startsWith("fixedStep")) {
                continue;
            }
            assertEquals(line1, line2);
        }

        reader1.close();
        reader2.close();


    }

    @Test
    public void testIndexedFasta() throws Exception {
        String fasta_file = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
        String infile = TestUtils.DATA_DIR + "bed/ecoli_out.test.bed";
        String outfile = TestUtils.DATA_DIR + "out/findextest.wig";
        String end_command = infile + " " + outfile + " " + fasta_file;
        String count_command = "count " + end_command;
        igvTools.run(count_command.split("\\s"));
        File outFi = new File(outfile);
        assertTrue(outFi.exists());
        outFi.delete();

        igvTools = new IgvTools();
        String index_command = "index " + end_command;
        igvTools.run(index_command.split("\\s"));
        String indexFile = infile + ".idx";
        File idxFi = new File(indexFile);
        assertTrue(idxFi.exists());
        idxFi.deleteOnExit();

        Genome genome = IgvTools.loadGenome(fasta_file, true);
        FeatureCodec codec = CodecFactory.getCodec(infile, genome);
        AbstractFeatureReader<Feature> reader = AbstractFeatureReader.getFeatureReader(infile, codec, true);
        String chr = "NC_000913_bb";
        Iterator<Feature> features = reader.query(chr, 5085, 5091);
        int count = 0;
        while (features.hasNext() && count < 100) {
            assertEquals(chr, features.next().getChr());
            count++;
        }
        assertEquals(3, count);

    }

    @Test
    public void testCountAllWF() throws Exception {
        String[] exclude = {"none", "density", "stddev", "count"};
        List<WindowFunction> wfList = new ArrayList<WindowFunction>();
        wfList.addAll(Arrays.asList(WindowFunction.values()));
        for (String s : exclude) {
            wfList.remove(WindowFunction.valueOf(s));
        }
        String inputFile = TestUtils.DATA_DIR + "bam/NA12878.SLX.sample.bam";
        tstCountWindowFunctions(inputFile, "All", wfList);
    }

    private void tstCountWindowFunctions(String inputFile, String chr, Iterable<WindowFunction> windowFunctions) throws Exception {
        String outputFile = TestUtils.DATA_DIR + "out/testCountWindowFunctions.tdf";


        String wfs = "";
        for (WindowFunction wf : windowFunctions) {
            wfs += wf.name() + ",";
        }
        wfs = wfs.substring(0, wfs.lastIndexOf(","));
        String[] cmd = {"count", "--windowFunctions", wfs, inputFile, outputFile, hg18id};
        igvTools.run(cmd);

        TDFReader reader = new TDFReader(new ResourceLocator(outputFile));
        for (WindowFunction wf : windowFunctions) {
            TDFDataset ds = reader.getDataset(chr, 0, wf);
            assertNotNull(ds);
        }
    }

    @Test
    //tile -z 1 --fileType mage-tab
    public void testTileMageTab() throws Exception {
        String mageTabFile = TestUtils.DATA_DIR + "mage-tab/test.data.txt";
        String outputFile = TestUtils.DATA_DIR + "mage-tab/test.data.tdf";
        String genfile = TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome";
        String command = "tile -z 1 --fileType mage-tab " + mageTabFile + " " + outputFile + " " + genfile;

        igvTools.run(command.split("\\s+"));


    }

}