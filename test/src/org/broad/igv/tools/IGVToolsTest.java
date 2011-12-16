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

import junit.framework.TestCase;
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;


public class IGVToolsTest extends TestCase {

    IgvTools igvTools;

    @Before
    public void setUp() throws Exception {
        super.setUp();
        igvTools = new IgvTools();
    }

    @After
    public void tearDown() throws Exception {
        igvTools = null;
        super.tearDown();
    }

    @Test
    public void testLinearIndex() throws IOException {

        String bedFile = "test/data/bed/test.bed";

        File idxFile = new File(bedFile + ".idx");
        if (idxFile.exists()) {
            idxFile.delete();
        }

        igvTools.doIndex(bedFile, 1, 16000);

        assertTrue(idxFile.exists());

        Index idx = IndexFactory.loadIndex(idxFile.getAbsolutePath());

        Block block = idx.getBlocks("chr1", 100, 200).get(0);
        assertEquals(0, block.getStartPosition());
        assertEquals(54, block.getSize());

        block = idx.getBlocks("chr1", 20000, 20001).get(0);
        assertEquals(54, block.getStartPosition());
        assertEquals(0, block.getSize());
    }

    @Test
    public void testIntervalIndex() throws IOException {

        // Create an interval tree index with 5 features per interval
        String testFile = "test/data/CEU.SRP000032.2010_03.genotypes.head.vcf";

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

        BasicFeatureSource bfr = BasicFeatureSource.getFeatureSource(testFile, new VCFCodec());
        Iterator<VariantContext> iter = bfr.query(chr, start, end);
        int count = 0;
        while (iter.hasNext()) {
            VariantContext feat = iter.next();
            int expStart = expectedStarts[count] - 1;
            assertEquals(expStart, feat.getStart());
            count++;
        }
        Assert.assertEquals(15, count);
    }


    @Test
    public void testLargeBed() throws Exception {
        //chr2:178,599,764-179,830,787 <- CONTAINS TTN
        String bedFile = "test/data/bed/Unigene.sorted.bed";

        // Interval index
        File indexFile = new File(bedFile + ".idx");
        if (indexFile.exists()) {
            indexFile.delete();
        }
        igvTools.doIndex(bedFile, IgvTools.INTERVAL_INDEX, 100);
        indexFile.deleteOnExit();


        String chr = "chr2";
        int start = 178599764;
        int end = 179830787;

        BasicFeatureSource bfr = BasicFeatureSource.getFeatureSource(bedFile, new BEDCodec());
        Iterator<BEDFeature> iter = bfr.query(chr, start, end);
        int count = 0;
        while (iter.hasNext()) {
            BEDFeature feat = iter.next();
            System.out.println(feat.getName());
            count++;
        }

    }


    @Test
    public void testVersion() throws IOException {
        String[] args = {"version"};
        //IgvTools.main(args);
        igvTools.run(args);
    }

    @Test
    public void testCompressWigFile() throws IOException {

        String inputFile = "test/data/phastCons_chr1.wig";
        testCompressOption(inputFile, 0, 0);
    }

    @Test
    public void testCompressCNFile() throws IOException {

        String inputFile = "test/data/HindForGISTIC.hg16.cn";
        testCompressOption(inputFile, 5000000, 6000000);
    }

    private void testCompressOption(String inputFile, int start, int end) throws IOException {
        String uncompressedFile = "uncompressed.tdf";
        String compressedFile = "compressed.tdf";
        String genome = "test/data/genomes/hg18.genome";


        String[] args = {"tile", "-z", "0", inputFile, uncompressedFile, genome};
        igvTools.run(args);

        args = new String[]{"tile", "-z", "0", "-c", inputFile, compressedFile, genome};
        (new IgvTools()).run(args);


        String dsName = "/chr1/raw";

        TDFDataset ds1 = TDFReader.getReader(uncompressedFile).getDataset(dsName);
        TDFDataset ds2 = TDFReader.getReader(compressedFile).getDataset(dsName);

        TDFTile t1 = ds1.getTiles(start, end).get(0);
        TDFTile t2 = ds2.getTiles(start, end).get(0);

        int nPts = t1.getSize();
        assertEquals(nPts, t2.getSize());

        for (int i = 0; i < nPts; i++) {
            assertEquals(t1.getStartPosition(i), t2.getStartPosition(i));
            assertEquals(t1.getValue(0, i), t2.getValue(0, i), 1.0e-6);
        }

        //(new File(uncompressedFile)).delete();
        //(new File(compressedFile)).delete();
    }

}