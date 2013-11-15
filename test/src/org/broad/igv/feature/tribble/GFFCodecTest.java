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
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureCodecParser;
import org.broad.igv.track.GFFFeatureSource;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Iterator;
import java.util.List;

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
}
