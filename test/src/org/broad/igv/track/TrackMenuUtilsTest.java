/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import com.google.common.base.Predicate;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
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
        ReferenceFrame.Range range = new ReferenceFrame.Range(chr, start, end);

        TrackMenuUtils.exportVisibleData(outPath, loadedTrack, range);

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
