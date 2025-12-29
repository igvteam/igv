package org.igv.feature.tribble;

import htsjdk.tribble.Feature;
import org.igv.AbstractHeadlessTest;
import org.igv.feature.BasicFeature;
import org.igv.feature.gff.GFFFeatureSource;
import org.igv.track.FeatureSource;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
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

        GFFFeatureSource src = new GFFFeatureSource(TribbleFeatureSource.getFeatureSource(new ResourceLocator(path), genome), GFFCodec.Version.GTF);

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
        FeatureSource source = new GFFFeatureSource(tribbleFeatureSource, GFFCodec.Version.GFF3);

        int featureCount = 0;
        Iterator<Feature> iter = source.getFeatures("chr7", 0, Integer.MAX_VALUE);
        while (iter.hasNext()) {
            BasicFeature bf = (BasicFeature) iter.next();
            featureCount++;
        }
        assertEquals(2, featureCount);

    }
}
