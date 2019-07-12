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

package org.broad.igv.util.collections;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureUtils;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 2/9/14
 *         Time: 9:38 PM
 */
public class CollUtilsTest {


    @Test
    public void testFilter() throws Exception {


        int nItems = 10000;
        List<Feature> features = new ArrayList(nItems);

        int start=0;
        int end=0;
        for (int i = 0; i < nItems; i++) {
            start = i * 5;
            end = start + 5;
            features.add(new BasicFeature("chr1", start, end));
        }

        List<Feature> filteredFeatures = CollUtils.filter(features, FeatureUtils.getOverlapPredicate("chr1", start, end));

        // We should find 1 feature (the last one)
        assertEquals(1, filteredFeatures.size());
        assertEquals(features.get(features.size()-1), filteredFeatures.get(0));
    }



}
