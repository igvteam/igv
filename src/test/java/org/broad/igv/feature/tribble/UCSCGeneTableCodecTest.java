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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/03/08
 */
public class UCSCGeneTableCodecTest extends AbstractHeadlessTest {

    @Test
    public void testLoadGenePred() throws Exception {
        String file = TestUtils.DATA_DIR + "gene/EnsembleGenes_sample.genepred";
        TrackLoader loader = new TrackLoader();
        ResourceLocator locator = new ResourceLocator(file);
        Genome genome = TestUtils.loadGenome();
        List<Track> tracks = loader.load(locator, genome);

        assertEquals(1, tracks.size());

        FeatureCodec codec = CodecFactory.getCodec(locator, genome);
        AbstractFeatureReader<Feature, ?> bfs = AbstractFeatureReader.getFeatureReader(locator.getPath(), codec, false);
        Iterable<Feature> iter = bfs.iterator();
        int count = 0;
        for (Feature f : iter) {
            if (count == 0) {
                assertEquals(67051161, f.getStart());
                assertEquals(67163158, f.getEnd());
            }
            count++;
        }

        assertEquals(74, count);

    }

    @Test
    public void testLoadUCSC() throws Exception {
        String file = TestUtils.DATA_DIR + "gene/UCSCgenes_sample.gene";
        TrackLoader loader = new TrackLoader();
        ResourceLocator locator = new ResourceLocator(file);
        Genome genome = TestUtils.loadGenome();
        List<Track> tracks = loader.load(locator, genome);

        assertEquals(1, tracks.size());

        FeatureCodec codec = CodecFactory.getCodec(locator.getPath(), genome);
        AbstractFeatureReader<Feature, ?> bfs = AbstractFeatureReader.getFeatureReader(locator.getPath(), codec, false);
        Iterable<Feature> iter = bfs.iterator();
        int count = 0;
        for (Feature f : iter) {
            if (count == 0) {
                assertEquals(1115, f.getStart());
                assertEquals(4121, f.getEnd());
            }
            count++;
        }

        assertEquals(382, count);

    }


}
