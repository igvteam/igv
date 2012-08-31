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

package org.broad.igv.variant;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.tribble.Feature;
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
