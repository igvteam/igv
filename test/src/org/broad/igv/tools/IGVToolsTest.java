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
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.junit.*;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;


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

        BasicFeatureSource bfr = BasicFeatureSource.getFeatureSource(testFile, codec);
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
    public void testCount() throws Exception {
        String inputFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        String outputFile = TestUtils.DATA_DIR + "/out/file_";
        String genome = TestUtils.DATA_DIR + "/genomes/hg18.unittest.genome";

        String[] opts = new String[]{"s", "sf", "sfb", "f", "b", ""};
        for (String opt : opts) {
            String[] args = {"count", "-a", "sc=" + opt, inputFile, outputFile + opt + ".tdf", genome};
            igvTools.run(args);

            TDFReader reader = TDFReader.getReader(outputFile + opt + ".tdf");
            assertTrue(reader.getDatasetNames().size() > 0);
        }
    }

}