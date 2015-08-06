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

package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.data.WiggleParser;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.*;
import java.util.*;

import static junit.framework.Assert.*;

public class IGVToolsCountTest extends AbstractHeadlessTest {

    IgvTools igvTools;

    private static final String hg18id = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";
    private static final int MAX_LINES_CHECK = 200;

    @Rule
    public TestRule testTimeout = new Timeout((int) 1e3 * 60 * 20);

    @Before
    public void setUp() throws Exception {
        super.setUp();
        igvTools = new IgvTools();

    }

    @After
    public void tearDown() throws Exception {
        super.tearDown();
        igvTools = null;
    }

    @Test
    public void testCountBEDoutTDF() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        tstCountStrandOpts(inputFile, "testtdf", "tdf", null, -1, -1);
    }

    @Test
    public void testCountBEDoutWig() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        tstCountStrandOpts(inputFile, "testwig", "wig", null, -1, -1);
        tstCountStrandMapOpts(inputFile, "testwig", "wig", null, -1, -1);
    }

    @Test
    public void testCountBEDoutWigCheckBW() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/test2.bed";
        String fullout = TestUtils.TMP_OUTPUT_DIR + "twig.wig";
        String input = "count " + inputFile + " " + fullout + " " + hg18id;
        String[] args = input.split("\\s+");
        igvTools.run(args);

        //Check file
        File outFile = new File(fullout);
        assertTrue(outFile.exists());
        assertTrue(outFile.canRead());

        WiggleParser parser = new WiggleParser(new ResourceLocator(fullout));
        WiggleDataset wgs = parser.parse();
        assertEquals(3, wgs.getChromosomes().length);
        for(String chr: new String[]{"chr1", "chr3", "chr7"}){
            assertNotNull(wgs.getData(null, chr));
        }

        BufferedReader reader = new BufferedReader(new FileReader(outFile));
        String line = reader.readLine();
        assertTrue(line.contains("wiggle_0"));
        //These shouldn't be present in the rest of the file
        while((line = reader.readLine()) != null){
            assertFalse(line.startsWith("#"));
            assertFalse(line.contains("wiggle_0"));
        }

    }

    @Test
    public void testCountBEDstdoutWig() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";

        //Write to a file
        String fileOutArg = TestUtils.DATA_DIR + "out/testfilewig.wig";
        String input = "count " + inputFile + " " + fileOutArg + " " + hg18id;
        String[] args = input.split("\\s+");
        igvTools.run(args);

        //Write to stdout. We redirect
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream oldOut = System.out;
        System.setOut(new PrintStream(os));

        input = "count " + inputFile + " " + IgvTools.STDOUT_FILE_STR + " " + hg18id;
        args = input.split("\\s+");
        (new IgvTools()).run(args);

        System.setOut(oldOut);

        BufferedReader fiReader = new BufferedReader(new FileReader(fileOutArg));
        BufferedReader strReader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(os.toByteArray())));

        String fiLine = null;
        String strLine = null;
        int lineCount = 0;
        while(true){
            fiLine = fiReader.readLine();
            strLine = strReader.readLine();

            assertEquals(fiLine, strLine);
            lineCount++;

            if(fiLine == null || strLine == null) break;
        }

        assertTrue(lineCount > 10);
    }

    @Test
    public void testCountBED_region() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        tstCountStrandOpts(inputFile, "testwig", "wig", "chr2", 178709699, 179008373);
        tstCountStrandMapOpts(inputFile, "testwig", "wig", "chr2", 178709699, 179008373);
    }

    @Test
    public void testCountSAM() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "sam/test_2.sam";
        tstCountStrandOpts(inputFile, "testwig", "wig", null, -1, -1);
    }

    @Test
    public void testCountBAM() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "bam/NA12878.SLX.sample.bam";
        tstCountStrandOpts(inputFile, "testwig", "wig", null, -1, -1);
        tstCountStrandMapOpts(inputFile, "testwig", "wig", null, -1, -1);
    }

    public void tstCountStrandOpts(String inputFile, String outputBase, String outputExt,
                                   String chr, int start, int end) throws Exception {

        String[] opts = new String[]{"--bases", "--strands=read", "--strands=first", "--strands=read --bases"};
        tstCountOptsGen(inputFile, outputBase, outputExt, chr, start, end, "", opts);
    }

    public void tstCountStrandMapOpts(String inputFile, String outputBase, String outputExt,
                                      String chr, int start, int end) throws Exception {

        String refOpt = "--minMapQuality 1";
        String[] opts = new String[]{"--strands=read " + refOpt, "--bases " + refOpt};
        tstCountOptsGen(inputFile, outputBase, outputExt, chr, start, end, refOpt, opts);
    }

    /**
     * Compare the output of count with the various options against the output of {@code refOpt}. That is,
     * we assume the row total for each will be equal.
     *
     * @param inputFile
     * @param outputBase
     * @param outputExt
     * @param chr
     * @param start
     * @param end
     * @param opts
     */
    public void tstCountOptsGen(String inputFile, String outputBase, String outputExt,
                                String chr, int start, int end, String refOpt, String[] opts) throws Exception {


        boolean query = chr != null && start >= 0 && end >= start + 1;
        String outputFile = TestUtils.TMP_OUTPUT_DIR + outputBase + "_";

        Map<String, float[]> rowTotals = new HashMap<String, float[]>(opts.length + 1);
        String[] allOpts = new String[opts.length + 1];
        allOpts[0] = refOpt;
        System.arraycopy(opts, 0, allOpts, 1, opts.length);

        for (int ind = 0; ind < allOpts.length; ind++) {
            String opt = allOpts[ind];
            int winsize = 5;
            if (query) {
                opt += " --windowSize " + winsize + " --query " + chr + ":" + start + "-" + end;
            }

            //In case we modify the option for querying as above
            if (ind == 0) refOpt = opt;

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
                WiggleDataset ds = (new WiggleParser(locator, IgvTools.loadGenome(hg18id))).parse();

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
            //Skip header lines
            if(line.startsWith("#") || line.contains("Step") || line.contains("wiggle_0")) continue;
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
    public void testCountDups() throws Exception {
        String inputFiname = "test_5duplicates";
        String ext = ".sam";
        String inputFile = TestUtils.DATA_DIR + "sam/" + inputFiname + ext;
        String outputFileND = TestUtils.TMP_OUTPUT_DIR + inputFiname + "_nodups" + ".tdf";
        String outputFileWithDup = TestUtils.TMP_OUTPUT_DIR + inputFiname + "_withdups" + ".tdf";

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
    public void testCountAliasedForward() throws Exception {
        String genfile = TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome";
        tstCountAliased(TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam",
                TestUtils.DATA_DIR + "sam/NA12878.muc1.test_modchr.sam",
                genfile);
    }

    @Test
    public void testCountAliasedReversed() throws Exception {
        String genfile = TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased_reversed.genome";
        tstCountAliased(TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam",
                TestUtils.DATA_DIR + "sam/NA12878.muc1.test_modchr.sam",
                genfile);
    }

    /**
     * Test count when using a custom alias file.
     * Should supply a file using normal chromosome names, and aliases
     * which are defined in the genome file.
     */
    public void tstCountAliased(String normfile, String aliasedfile, String genfile) throws Exception {

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
    public void testCountIndexedFasta() throws Exception {
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

        Genome genome = IgvTools.loadGenome(fasta_file);
        FeatureCodec codec = CodecFactory.getCodec(infile, genome);
        AbstractFeatureReader<Feature, ?> reader = AbstractFeatureReader.getFeatureReader(infile, codec, true);
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

}