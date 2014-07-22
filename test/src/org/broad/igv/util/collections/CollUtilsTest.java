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
