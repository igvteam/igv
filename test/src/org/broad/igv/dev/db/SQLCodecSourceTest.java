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

package org.broad.igv.dev.db;

import junit.framework.Assert;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;

public class SQLCodecSourceTest {

    @Test
    public void testLoadBED() throws Exception {

        AsciiFeatureCodec codec = new IGVBEDCodec();

        String host = (new File(TestUtils.DATA_DIR)).getAbsolutePath();
        String path = "sql/unigene.db";
        String url = DBManager.createConnectionURL("sqlite", host, path, null);
        ResourceLocator locator = new ResourceLocator(url);
        String table = "unigene";


        SQLCodecSource reader = new SQLCodecSource(locator, codec, table, "chrom", "chromStart", "chromEnd", 1, Integer.MAX_VALUE);
        Iterator<Feature> SQLFeatures = reader.iterator();

        String bedFile = host + "/bed/Unigene.sample.bed";
        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec, false);
        Iterator<Feature> fileFeatures = bfr.iterator();

        int count = 0;
        while (SQLFeatures.hasNext()) {
            Feature f = SQLFeatures.next();
            Feature fileFeature = fileFeatures.next();
            Assert.assertEquals(fileFeature.getChr(), f.getChr());
            Assert.assertEquals(fileFeature.getStart(), f.getStart());
            Assert.assertEquals(fileFeature.getEnd(), f.getEnd());
            count++;
        }

        Assert.assertEquals(72, count);
    }
}