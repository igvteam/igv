/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.tribble;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
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

        FeatureCodec codec = CodecFactory.getCodec(locator.getPath(), genome);
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(locator.getPath(), codec, false);
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
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(locator.getPath(), codec, false);
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
