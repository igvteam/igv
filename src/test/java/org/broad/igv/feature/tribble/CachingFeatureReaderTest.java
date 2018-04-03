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

package org.broad.igv.feature.tribble;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 *
 */
public class CachingFeatureReaderTest {

    String path = TestUtils.LARGE_DATA_DIR + "CEU.SRP000032.2010_03.genotypes.vcf.gz";
    AbstractFeatureReader baseReader;
    CachingFeatureReader cacheReader;

    @Before
    public void setUp() throws IOException {
        Genome genome = null; // <= don't do chromosome alias conversion
        FeatureCodec codec = CodecFactory.getCodec(path, null);
        baseReader =    AbstractFeatureReader.getFeatureReader(path, path + ".tbi", codec, true);
        cacheReader = new CachingFeatureReader(baseReader);

    }

    @Test
    /*
    1 182773022 182773023
1 182773812 182773813
1 182774171 182774172
1 182774425 182774426
1 182774496 182774497
1 182774630 182774631
1 182774894 182774895
1 182775195 182775196
1 182775253 182775254
1 182776704 182776705
1 182776741 182776742
     */
    public void testGetSequenceNames() throws Exception {

        List<String> seqNames = cacheReader.getSequenceNames();
        assertEquals(23, seqNames.size());
    }

    @Test
    public void testQuery() throws Exception {

        Set<String> baseReaderLoci = new HashSet();
        final int start = 182773022;
        final int end = 182776742;
        final String chr = "1";
        Iterator<VCFVariant> iter = baseReader.query(chr, start, end);

        while (iter.hasNext()) {
            VCFVariant var = iter.next();
            assertEquals(chr, var.getChr());
            assertTrue(var.getEnd() >= start && var.getStart() <= end);
            baseReaderLoci.add(getLocusString(var));
        }
        assertTrue(baseReaderLoci.size() > 0);

        // Now use CachingFeatureReader and insure results are the same
        Set<String> cacheReaderLoci = new HashSet();
        iter = cacheReader.query(chr, start, end);
        while (iter.hasNext()) {
            VCFVariant var = iter.next();
            assertEquals(chr, var.getChr());
            assertTrue(var.getEnd() >= start && var.getStart() <= end);
            cacheReaderLoci.add(getLocusString(var));
        }

        assertEquals(baseReaderLoci.size(), cacheReaderLoci.size());
        for (String locus : cacheReaderLoci) {
            assertTrue(baseReaderLoci.contains(locus));
        }
    }

    String getLocusString(VCFVariant variant) {
        return variant.getChr() + ":" + variant.getStart() + "-" + variant.getEnd();
    }
}
