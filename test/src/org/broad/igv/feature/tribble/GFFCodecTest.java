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
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.GFFFeatureSource;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2013-Mar-21
 */
public class GFFCodecTest extends AbstractHeadlessTest {

    /**
     * Make sure we parse the attributes to get the name of this feature
     * GTF has a bunch of different ones
     *
     * @throws Exception
     */
    @Test
    public void testGetNameGTF() throws Exception {
        String path = TestUtils.DATA_DIR + "gtf/transcript_id.gtf";
        String expName = "YAL069W";

        GFFFeatureSource src = new GFFFeatureSource(TribbleFeatureSource.getFeatureSource(new ResourceLocator(path), null));

        Iterator<Feature> iter = src.getFeatures("I", 0, Integer.MAX_VALUE);
        while (iter.hasNext()) {
            BasicFeature bf = (BasicFeature) iter.next();
            assertEquals(expName, bf.getName());
        }

    }

    /**
     * Insure we can parse a GFF file that includes a fasta section.
     *
     * @throws Exception
     */
    @Test
    public void testGFFWithFasta() throws Exception {

        String path = TestUtils.DATA_DIR + "gff/gffWithFasta.gff";
        final ResourceLocator locator = new ResourceLocator(path);

        TribbleFeatureSource tribbleFeatureSource = TribbleFeatureSource.getFeatureSource(locator, genome);
        FeatureSource source = new GFFFeatureSource(tribbleFeatureSource);

        int featureCount = 0;
        Iterator<Feature> iter = source.getFeatures("chr7", 0, Integer.MAX_VALUE);
        while (iter.hasNext()) {
            BasicFeature bf = (BasicFeature) iter.next();
            featureCount++;
        }
        assertEquals(2, featureCount);

    }
}
