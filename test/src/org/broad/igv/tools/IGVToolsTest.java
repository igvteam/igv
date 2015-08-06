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
import org.broad.igv.data.Dataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.FastaIndex;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.sam.reader.SamUtils;
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.igv.tools.sort.SorterTest;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
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

    private String doStandardIndex(String inputFile, String expectedExtension) throws IOException {
        String indDir = TestUtils.TMP_OUTPUT_DIR;
        TestUtils.clearOutputDir();

        String indPath = igvTools.doIndex(inputFile, indDir, IgvTools.LINEAR_INDEX, IgvTools.LINEAR_BIN_SIZE);
        File indFile = new File(indPath);

        //Check that only the index file we intended exists
        assertTrue(indFile.exists());
        assertTrue(indPath.endsWith(expectedExtension));

        final Set<String> exts = new HashSet<String>();
        for (String ext : new String[]{".idx", ".sai", ".bai", ".fai"}) {
            exts.add(ext);
        }
        File[] files = (new File(indDir)).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return exts.contains((Preprocessor.getExtension(name)));
            }
        });
        assertEquals("Extra files in output directory", 1, files.length);

        return indFile.getAbsolutePath();
    }


    @Test
    public void testIndexSam() throws Exception {
        String samFile = TestUtils.DATA_DIR + "sam/NA12878.muc1.test2.sam";
        String samFileIdx = doStandardIndex(samFile, "sai");

        FeatureIndex idx = SamUtils.getIndexFor(samFile);
        assertTrue(idx.containsChromosome("chr1"));
        assertEquals(1, idx.getIndexedChromosomes().size());
    }

    @Test
    public void testIndexFasta() throws Exception {
        String inFile = TestUtils.DATA_DIR + "fasta/ecoli_out.padded2.fasta";
        String indPath = doStandardIndex(inFile, "fai");

        FastaIndex index = new FastaIndex(indPath);
        assertEquals(1, index.getSequenceNames().size());
        assertNotNull(index.getIndexEntry("NC_000913_bb"));
    }

    @Test
    public void testLinearIndex() throws IOException {

        String bedFile = TestUtils.DATA_DIR + "bed/test.bed";

        String idxPath = doStandardIndex(bedFile, "idx");

        Index idx = IndexFactory.loadIndex(idxPath);

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
        igvTools.doIndex(testFile, null, 2, 5);
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
        String inputFile = TestUtils.DATA_DIR + "wig/phastCons_chr1_small.wig";
        testTile(inputFile, 300, 1300);
    }

    @Test
    public void testTileCNFile() throws IOException {
        String inputFile = TestUtils.DATA_DIR + "cn/HindForGISTIC.hg16.cn";
        testTile(inputFile, 5000000, 5500000);
    }


    @Test
    public void testTileGCT_01() throws IOException {
        String inputFile = TestUtils.DATA_DIR + "gct/OV.transcriptome__agilentg4502.data.txt";
        String outFilePath = TestUtils.DATA_DIR + "out/testTileGCT.wig";
        String[] args = {"tile", "-z", "1", "--fileType", "mage-tab", inputFile, outFilePath, hg18id};
        igvTools.run(args);
    }

    @Test
    public void testTileGCT_02() throws IOException {
        String inputFile = TestUtils.DATA_DIR + "gct/GBM.methylation__sampled.data.txt";
        String outFilePath = TestUtils.DATA_DIR + "out/testTileGCT.wig";
        String[] args = new String[]{"tile", "-z", "1", "--fileType", "mage-tab", inputFile, outFilePath, hg18id};
        igvTools.run(args);

    }


    private void testTile(String inputFile, int start, int end) throws IOException {
        String file1 = TestUtils.DATA_DIR + "out/file1.tdf";
        String file2 = TestUtils.DATA_DIR + "out/file2.tdf";

        //todo Compare 2 outputs more meaningfully
        String[] args = {"toTDF", "-z", "1", "--windowFunctions", "min", inputFile, file1, hg18id};
        igvTools.run(args);

        FeatureDB.clearFeatures();
        Runtime.getRuntime().gc();

        args = new String[]{"toTDF", "-z", "1", "--windowFunctions", "max", inputFile, file2, hg18id};
        igvTools.run(args);


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

    /**
     * Test iterating through a merged bam file(actually a list with the same file duplicated).  Each record
     * should appear twice (since the list contains the same file twice), and in coordinate sort order.
     *
     * @throws Exception
     */
    @Test
    public void testIterateMergedBam() throws Exception {
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

    @Test
    public void testSort() throws Exception {
        tstSort(false);
    }

    @Test
    public void testSortStdout() throws Exception {
        tstSort(true);
    }

    public void tstSort(boolean writeToStdOut) throws Exception {
        String inputFiname = "Unigene.unsorted.bed";
        String inputFile = TestUtils.DATA_DIR + "bed/" + inputFiname;
        String outputFile = TestUtils.TMP_OUTPUT_DIR + inputFiname + ".sorted";
        File oFile = new File(outputFile);
        oFile.deleteOnExit();
        String outputArg = outputFile;

        //This looks a bit funny, but for ease of testing we redirect stdout to a file
        //Mostly just concerned about spurious log statements getting into the file
        if(writeToStdOut){
            System.setOut(new PrintStream(new FileOutputStream(oFile)));
            outputArg = IgvTools.STDOUT_FILE_STR;
        }

        String input = "sort --tmpDir=./ --maxRecords=50 " + inputFile + " " + outputArg;
        igvTools.run(input.split("\\s+"));

        int numlines = SorterTest.checkFileSorted(oFile, 0, 1, 0);
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
        String outputFile = TestUtils.TMP_OUTPUT_DIR + inputFiname + "_formatted" + ext;
        File oFile = new File(outputFile);
        oFile.deleteOnExit();

        String input = "formatexp " + inputFile + " " + outputFile;
        igvTools.run(input.split("\\s+"));
        Genome genome = TestUtils.loadGenome();

        ExpressionFileParser parser = new ExpressionFileParser(new ResourceLocator(outputFile), null, genome);
        Dataset ds = parser.createDataset();
        assertEquals(10, ds.getChromosomes().length);
    }

    @Ignore("Missing data file")
    @Test
    public void testTileMageTab() throws Exception {
        String mageTabFile = TestUtils.DATA_DIR + "mage-tab/test.data.txt";
        String outputFile = TestUtils.DATA_DIR + "mage-tab/test.data.tdf";
        String genfile = TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome";
        String command = "tile -z 1 --fileType mage-tab " + mageTabFile + " " + outputFile + " " + genfile;

        igvTools.run(command.split("\\s+"));
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

}