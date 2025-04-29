package org.broad.igv.feature;

import htsjdk.tribble.Feature;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class DataURLParserTest extends AbstractHeadlessTest {


    @Test
    public void testDataURL() throws IOException {

        String f = TestUtils.DATA_DIR + "bed/H3K4me1_sample_bed6.bed";
        String dataURL = DataURLParser.createDataURL(new File(f));

        DataURLParser parser = new DataURLParser();

        parser.parseFeatures(dataURL, "bed", genome);
        List<Feature> features = parser.getFeatures();
        assertEquals(18, features.size());

        TrackProperties trackProperties = parser.getTrackProperties();
        assertEquals("H3k4me1", trackProperties.getName());

        // chr1	226823	226846	U0	0	+
        BasicFeature lastFeature = (BasicFeature) features.get(17);
        assertEquals("chr1", lastFeature.getChr());
        assertEquals(226823, lastFeature.getStart());
        assertEquals(Strand.POSITIVE, lastFeature.getStrand());

    }


}