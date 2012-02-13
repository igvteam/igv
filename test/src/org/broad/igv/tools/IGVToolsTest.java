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

package org.broad.igv.tools;

import org.broad.igv.Globals;
import org.broad.igv.data.Dataset;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.data.WiggleParser;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.sam.reader.SamUtils;
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.igv.tools.sort.SorterTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.junit.*;

import java.io.*;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.*;


public class IGVToolsTest {

    IgvTools igvTools;

    @Before
    public void setUp() throws Exception {
        Globals.setHeadless(true);
        igvTools = new IgvTools();
    }

    @After
    public void tearDown() throws Exception {
        igvTools = null;
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        //TestUtils.clearOutputDir();
    }

    @Test
    public void testIndexSam() throws Exception {
        String samFile = TestUtils.DATA_DIR + "/sam/NA12878.muc1.test2.sam";
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

        String bedFile = TestUtils.DATA_DIR + "/bed/test.bed";

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
        assertEquals("Unexpected block size", 100, block.getSize());

        List<Block> allblocks = idx.getBlocks("chr1", 1, Integer.MAX_VALUE);
        //5 lines, get broken up into 2 blocks
        assertEquals(2, allblocks.size());
    }

    @Test
    public void testIntervalIndex33() throws Exception {
        String testFile = TestUtils.LARGE_DATA_DIR + "/CEU.SRP000032.2010_03_v3.3.genotypes.head.vcf";
        FeatureCodec codec = new VCF3Codec();
        testIntervalIndex(testFile, codec);
    }

    @Test
    public void testIntervalIndex40() throws Exception {
        String testFile = TestUtils.LARGE_DATA_DIR + "/CEU.SRP000032.2010_03_v4.0.genotypes.head.vcf";
        FeatureCodec codec = new VCFCodec();
        testIntervalIndex(testFile, codec);
    }

    public void testIntervalIndex(String testFile, FeatureCodec codec) throws IOException {

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

        String inputFile = TestUtils.LARGE_DATA_DIR + "/phastCons_chr1.wig";
        testTile(inputFile, 0, 0);
    }

    @Test
    public void testTileCNFile() throws IOException {

        String inputFile = TestUtils.DATA_DIR + "/cn/HindForGISTIC.hg16.cn";
        testTile(inputFile, 5000000, 6000000);
    }

    private void testTile(String inputFile, int start, int end) throws IOException {
        String file1 = TestUtils.DATA_DIR + "/out/file1.tdf";
        String file2 = TestUtils.DATA_DIR + "/out/file2.tdf";
        String genome = TestUtils.DATA_DIR + "/genomes/hg18.unittest.genome";


        //todo Compare 2 outputs more meaningfully
        String[] args = {"tile", "-z", "1", "--windowFunctions", "min", inputFile, file1, genome};
        igvTools.run(args);

        args = new String[]{"tile", "-z", "2", "--windowFunctions", "max", inputFile, file2, genome};
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
        tstCount("testtdf", "tdf", null, -1, -1);
    }

    @Test
    public void testCountWIG() throws Exception {
        tstCount("testwig", "wig", null, -1, -1);
        tstCount("testwig", "wig", "chr2", 178709699, 179008373);
    }

    public void tstCount(String outputBase, String outputExt, String chr, int start, int end) throws Exception {
        String inputFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        String outputFile = TestUtils.DATA_DIR + "/out/" + outputBase + "_";
        String genome = TestUtils.DATA_DIR + "/genomes/hg18.unittest.genome";

        boolean query = chr != null && start >= 0 && end >= start + 1;

        String[] opts = new String[]{"", "--bases", "--strand=read", "--strand=first", "--strand=read --bases"};
        for (int ind = 0; ind < opts.length; ind++) {
            String opt = opts[ind];
            int winsize = 5;
            if (query) {
                opt += " --windowSize " + winsize + " --query " + chr + ":" + start + "-" + end;
            }

            String fullout = outputFile + ind + "." + outputExt;
            String input = "count " + opt + " " + inputFile + " " + fullout + " " + genome;
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
                        }
                    }
                }
            } else {
                File outFile = new File(fullout);
                assertTrue(outFile.exists());
                assertTrue(outFile.canRead());

                ResourceLocator locator = new ResourceLocator(fullout);
                WiggleDataset ds = (new WiggleParser(locator, IgvTools.loadGenome(genome, true))).parse();

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
    }

    @Test
    public void testCountBAMList() throws Exception {
        String listPath = TestUtils.DATA_DIR + "/bam/2largebams.bam.list";
        File listFile = new File(listPath);
        listFile.delete();
        listFile.deleteOnExit();
        //We generate the file on each test, because largedata dir can change
        String[] largebams = new String[]{"HG00171.hg18.bam", "HG00171.hg18.bam"};
        for (int ii = 0; ii < largebams.length; ii++) {
            largebams[ii] = TestUtils.LARGE_DATA_DIR + "/" + largebams[ii];
        }
        FileWriter writer = new FileWriter(listFile);
        for (String s : largebams) {
            writer.write(s + "\n");
        }
        writer.close();

        //Test list file
        tstCountBamList(listPath);

        //Test comma-separated list
        String listArg = largebams[0];
        for (int ss = 1; ss < largebams.length; ss++) {
            listArg += "," + largebams[ss];
        }
        tstCountBamList(listArg);

    }

    private void tstCountBamList(String listArg) throws Exception {
        String outputFile = TestUtils.DATA_DIR + "/out/file_";
        String genome = TestUtils.DATA_DIR + "/genomes/hg18.unittest.genome";

        String[] opts = new String[]{"--strand=read", "--strand=first", ""};

        for (int ind = 0; ind < opts.length; ind++) {
            String opt = opts[ind];
            String fullout = outputFile + ind + ".tdf";
            String input = "count " + opt + " " + listArg + " " + fullout + " " + genome;
            String[] args = input.split("\\s+");
            igvTools.run(args);

            TDFReader reader = TDFReader.getReader(fullout);
            assertTrue(reader.getDatasetNames().size() > 0);
        }
    }

    @Test
    public void testSort() throws Exception {
        String inputFiname = "Unigene.unsorted.bed";
        String inputFile = TestUtils.DATA_DIR + "/bed/" + inputFiname;
        String outputFile = TestUtils.DATA_DIR + "/out/" + inputFiname + ".sorted";
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
        String inputFile = TestUtils.DATA_DIR + "/gct/" + inputFiname + ext;
        String outputFile = TestUtils.DATA_DIR + "/out/" + inputFiname + "_formatted" + ext;
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
    public void testCountAliased() throws Exception {
        tstCountAliased(TestUtils.DATA_DIR + "/sam/NA12878.muc1.test.sam",
                TestUtils.DATA_DIR + "/sam/NA12878.muc1.test_modchr.sam");
    }


    /**
     * Test count when using a custom alias file.
     * Should supply a file using normal chromosome names, and aliases
     * which are defined in the genome file. These files are expected to
     * match exactly
     */
    public void tstCountAliased(String normfile, String aliasedfile) throws Exception {
        String genfile = TestUtils.DATA_DIR + "/genomes/hg18_truncated_aliased.genome";
        String outfile = TestUtils.DATA_DIR + "/out/tmpcount1.wig";
        File outFi = new File(outfile);
        outFi.delete();
        //Count aliased file
        String command = "count " + normfile + " " + outfile + " " + genfile;
        igvTools.run(command.split("\\s+"));

        //Count non-aliased file
        String outfile2 = TestUtils.DATA_DIR + "/out/tmpcount2.wig";
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
            assertEquals(line1, line2);
        }


    }

    @Test
    public void testIndexedFasta() throws Exception {
        String fasta_file = TestUtils.DATA_DIR + "/fasta/ecoli_out.padded.fasta";
        String infile = TestUtils.DATA_DIR + "/bed/ecoli_out.test.bed";
        String outfile = TestUtils.DATA_DIR + "/out/findextest.wig";
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
    }

}