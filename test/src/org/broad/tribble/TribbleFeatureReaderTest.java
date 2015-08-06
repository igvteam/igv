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

package htsjdk.tribble;


import org.broad.igv.util.TestUtils;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * User: jrobinso
 * Date: Jul 17, 2010
 * Time: 9:28:49 PM
 */
public class TribbleFeatureReaderTest {

    static String testFile = TestUtils.LARGE_DATA_DIR + "CEU.SRP000032.2010_03_v4.0.genotypes.head.vcf";
    static AbstractFeatureReader<VariantContext, ?> bfr;

    @BeforeClass
    public static void setUpClass() throws IOException {
        File idxFile = new File(testFile + ".idx");
        FeatureCodec codec = new VCFCodec();
        createIndex(idxFile, codec);

        bfr = AbstractFeatureReader.getFeatureReader(testFile, codec);
    }

    @AfterClass
    public static void tearDownClass() throws IOException {
        bfr.close();
    }

    @Test
    public void testQuery() throws IOException {

        String chr = "1";
        int start = 1718546;
        int end = 1748915;
        int[] expectedStarts = {1718547, 1718829, 1723079, 1724830, 1731376, 1733967, 1735586, 1736016, 1738594,
                1739272, 1741124, 1742815, 1743224, 1748886, 1748914};


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
    public void testGetSequenceNames() throws Exception {
        Set<String> expectedSequences = new HashSet(Arrays.asList("1", "2"));

        int count = 0;
        for (String s : bfr.getSequenceNames()) {
            assertTrue(expectedSequences.contains(s));
            count++;
        }

        assertEquals(expectedSequences.size(), count);
    }


    private static void createIndex(File idxFile, FeatureCodec codec) throws IOException {
        if (idxFile.exists()) {
            idxFile.delete();
        }

        // Create the index  
        Index idx = IndexFactory.createIntervalIndex(new File(testFile), codec, 10);

        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        } finally {
            if (stream != null) {
                stream.close();
            }
        }
        idxFile.deleteOnExit();

    }

}