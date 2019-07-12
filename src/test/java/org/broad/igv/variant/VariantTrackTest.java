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

package org.broad.igv.variant;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.tribble.Feature;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Aug-29
 */
public class VariantTrackTest extends AbstractHeadlessTest {

    @Test
    public void testLoadMultiAlleleFreqs() throws Exception {
        String filePath = TestUtils.DATA_DIR + "vcf/multi_allele_freqs.vcf";
        TestUtils.createIndex(filePath);

        Track newTrack = (new TrackLoader()).load(new ResourceLocator(filePath), genome).get(0);
        VariantTrack variantTrack = (VariantTrack) newTrack;

        List<Feature> featuresList = variantTrack.getFeatures("chr1", 542939 - 10, 543702 + 10);
        assertEquals(6, featuresList.size());

        double[][] expFreqs = {{0.067}, {0.067}, {0.143, 0.786}, {0.067}, {0.067}, {-1.0}};
        assert featuresList.size() == expFreqs.length;
        for (int fi = 0; fi < expFreqs.length; fi++) {
            double[] expFreq = expFreqs[fi];
            double[] actFreq = ((VCFVariant) featuresList.get(fi)).getAlleleFreqs();
            Assert.assertArrayEquals(expFreq, actFreq, 1e-12);
        }
    }

}
