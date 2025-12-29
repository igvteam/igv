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
