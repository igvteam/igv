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

package org.broad.igv.track;

import com.google.common.base.Predicate;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Sep-26
 */
public class TrackMenuUtilsTest extends AbstractHeadlessTest {

    @Test
    public void testExportData_BED() throws Exception{
        TrackLoader loader = new TrackLoader();
        ResourceLocator locator = new ResourceLocator(TestUtils.DATA_DIR + "bed/test.bed");
        List<Track> loadedTrack = loader.load(locator, genome);

        String outPath = TestUtils.TMP_OUTPUT_DIR + "tmpOut.bed";

        String chr = "chr1";
        int start = 150;
        int end = 350;

        Predicate<Feature> overlapPred = FeatureUtils.getOverlapPredicate(chr, start, end);
        Range range = new Range(chr, start, end);

        TrackMenuUtils.exportVisibleFeatures(outPath, loadedTrack, range);

        AbstractFeatureReader bfs = AbstractFeatureReader.getFeatureReader(outPath, CodecFactory.getCodec(outPath, genome), false);
        Iterator<Feature> iter = bfs.iterator();
        int count = 0;
        while(iter.hasNext()){
            Feature feat = iter.next();

            assertTrue(overlapPred.apply(feat));

            count += 1;
        }
        assertEquals(2, count);
    }
}
