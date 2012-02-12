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

package org.broad.tribble;


import org.broad.igv.util.TestUtils;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
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

    static String testFile = TestUtils.LARGE_DATA_DIR + "/CEU.SRP000032.2010_03_v4.0.genotypes.head.vcf";
    static AbstractFeatureReader<VariantContext> bfr;

    @BeforeClass
    public static void setUp() throws IOException {
        File idxFile = new File(testFile + ".idx");
        FeatureCodec codec = new VCFCodec();
        createIndex(idxFile, codec);

        bfr = AbstractFeatureReader.getFeatureReader(testFile, codec);
    }

    @AfterClass
    public static void tearDown() throws IOException {
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
            int expStart = expectedStarts[count];// - 1;
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

    }
}