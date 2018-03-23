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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Aug 9, 2010
 */
public class TribbleIndexTest extends AbstractHeadlessTest {

    /**
     * chr2	1	200000000	LONG_FEATURE
     * ...
     * chr2	179098961	179380395	Hs.134602
     * chr2	179209546	179287210	Hs.620337
     * chr2	179266309	179266748	Hs.609465
     * chr2	179296428	179300012	Hs.623987
     * chr2	179302952	179303488	Hs.594545
     */

    @Test
    public void testLinearIndex() throws Exception {

        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179266309 - 1;
        int end = 179303488 + 1;
        int expectedCount = 6;

        Set<String> expectedNames = new HashSet<String>(Arrays.asList("Hs.134602", "Hs.620337", "Hs.609465", "Hs.623987",
                "Hs.594545", "LONG_FEATURE"));

        // Interval index
        TestUtils.createIndex(bedFile, IgvTools.LINEAR_INDEX, 500);

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, new IGVBEDCodec());
        Iterator<BasicFeature> iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BasicFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(expectedNames.contains(feature.getName()));
            countInterval++;
        }

        assertEquals(expectedCount, countInterval);

    }

    @Test
    /**
     * Test interval tree index
     * chr2	179098961	179380395	Hs.134602
     * chr2	179209546	179287210	Hs.620337
     * chr2	179266309	179266748	Hs.609465
     * chr2	179296428	179300012	Hs.623987
     * chr2	179302952	179303488	Hs.594545
     *
     */
    public void testIntervalTree() throws Exception {
        //chr2:179,222,066-179,262,059<- CONTAINS TTN
        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179266309 - 1;
        int end = 179303488 + 1;
        int expectedCount = 6;

        Set<String> expectedNames = new HashSet<String>(Arrays.asList("Hs.134602", "Hs.620337", "Hs.609465", "Hs.623987",
                "Hs.594545", "LONG_FEATURE"));

        // Interval index
        TestUtils.createIndex(bedFile, IgvTools.INTERVAL_INDEX, 1);

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, new IGVBEDCodec());
        Iterator<BasicFeature> iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BasicFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(expectedNames.contains(feature.getName()));
            countInterval++;
        }

        assertEquals(expectedCount, countInterval);
    }

    @Test
    public void testReadSingleVCF() throws Exception {
        String file = TestUtils.DATA_DIR + "vcf/indel_variants_onerow.vcf";
        String chr = "chr9";
        // Linear index
        TestUtils.createIndex(file);

        // First test query
        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(file, new VCFCodec());
        Iterator<htsjdk.variant.variantcontext.VariantContext> iter = bfr.query(chr, 5073767 - 5, 5073767 + 5);
        int count = 0;
        while (iter.hasNext()) {
            htsjdk.variant.variantcontext.VariantContext feat = iter.next();
            assertEquals("chr9", feat.getChr());
            assertEquals(feat.getStart(), 5073767);
            assertTrue(feat.hasAttribute("MapQs"));
            count++;
        }
        assertEquals(1, count);

        // Test non-indexed access (iterator)
        iter = bfr.iterator();
        count = 0;
        while (iter.hasNext()) {
            htsjdk.variant.variantcontext.VariantContext feat = iter.next();
            assertEquals("chr9", feat.getChr());
            assertEquals(feat.getStart(), 5073767);
            assertTrue(feat.hasAttribute("MapQs"));
            count++;
        }
        assertEquals(1, count);

        //Do similar as above, but have a different test file
        file = TestUtils.DATA_DIR + "vcf/outputPileup.flt1.vcf";
        chr = "1";
        // Linear index
        TestUtils.createIndex(file);

        bfr = AbstractFeatureReader.getFeatureReader(file, new VCFCodec());
        iter = bfr.query(chr, 984163 - 5, 984163 + 5);
        count = 0;
        while (iter.hasNext()) {
            htsjdk.variant.variantcontext.VariantContext feat = iter.next();
            assertEquals(chr, feat.getChr());
            if (count == 0) {
                assertEquals(984163, feat.getStart());
            }
            count++;
        }
        assertEquals(1, count);

    }


}
