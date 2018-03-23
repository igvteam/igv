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

package org.broad.igv.dev.db;

import junit.framework.Assert;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

public class SQLCodecSourceTest {

    private SQLCodecSource getUnigene(String path) {
        AsciiFeatureCodec codec = new IGVBEDCodec();
        String host = (new File(TestUtils.DATA_DIR)).getAbsolutePath();

        String url = DBManager.createConnectionURL("sqlite", host, path, null);
        ResourceLocator locator = new ResourceLocator(url);
        String tableName = "unigene";

        DBProfile.DBTable table = new DBProfile.DBTable(locator, tableName, "n/a", null, "chrom", "chromStart", "chromEnd", 1, Integer.MAX_VALUE, null, null, null);
        SQLCodecSource reader = new SQLCodecSource(table, codec);
        return reader;
    }

    //Check that querying returns sorted features
    @Test
    public void testQueryBEDUnsorted() throws Exception {
        String path = "sql/Unigene.unsorted.db";
        SQLCodecSource reader = getUnigene(path);
        reader.setFeatureWindowSize(Integer.MAX_VALUE / 2);
        Iterator<Feature> features = reader.getFeatures("chr2", 0, Integer.MAX_VALUE / 4);
        assertNotNull(features);
        int count = TestUtils.assertFeatureIteratorSorted(features);
        assertEquals(71, count);
    }

    //Check that iterating returns sorted features
    @Test
    public void testIterateBEDUnsorted() throws Exception {
        String path = "sql/Unigene.unsorted.db";
        SQLCodecSource reader = getUnigene(path);
        Iterator<Feature> features = reader.iterator();
        int count = TestUtils.assertFeatureIteratorSorted(features);
        assertEquals(71, count);
    }

    @Test
    public void testLoadBED() throws Exception {

        AsciiFeatureCodec codec = new IGVBEDCodec();
        String path = "sql/unigene.db";

        SQLCodecSource reader = getUnigene(path);
        Iterator<Feature> SQLFeatures = reader.iterator();

        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec, false);
        Iterator<Feature> fileFeatures = bfr.iterator();

        int count = 0;
        while (SQLFeatures.hasNext()) {
            Feature f = SQLFeatures.next();
            Feature fileFeature = fileFeatures.next();
            TestUtils.assertFeaturesEqual(fileFeature, f);
            count++;
        }

        Assert.assertEquals(72, count);
    }

    //Don't support reordering by index
    //@Test
    public void testLoadReorderedColumnsIndex() throws Exception {
        String profilePath = TestUtils.DATA_DIR + "sql/unsorted.colsreordered.byindex.dbxml";
        tstLoadReorderedColumns(profilePath);
    }

    @Test
    public void testLoadReorderedColumnsLabel() throws Exception {
        String profilePath = TestUtils.DATA_DIR + "sql/unsorted.colsreordered.bylabel.dbxml";
        tstLoadReorderedColumns(profilePath);
    }

    public void tstLoadReorderedColumns(String profilePath) throws Exception {
        SQLCodecSource source0 = SQLCodecSource.getFromProfile(profilePath, "unigene");
        assertNotNull(source0);

        SQLCodecSource source1 = getUnigene("sql/Unigene.unsorted.db");

        Iterator<Feature> features0 = source0.iterator();
        Iterator<Feature> features1 = source1.iterator();

        while (features0.hasNext()) {
            Feature act_feat = features0.next();
            Feature exp_feat = features1.next();
            TestUtils.assertFeaturesEqual(exp_feat, act_feat);
        }

    }


    private SQLCodecSource getVCFsrc() throws Exception {
        String profilePath = TestUtils.DATA_DIR + "sql/SRP32_v40.dbxml";
        return SQLCodecSource.getFromProfile(profilePath, "SRP32_v40");

    }

    @Test
    public void testIterateVCF() throws Exception {

        int featCount = 0;
        Iterator<Feature> iterator = getVCFsrc().iterator();
        while (iterator.hasNext()) {
            iterator.next();
            featCount++;
        }

        assertEquals(3963, featCount);
    }

    @Test
    public void testQueryVCF() throws Exception {
        SQLCodecSource source = getVCFsrc();
        //Target row:
        //1	201813161	rs2248941	T	C	.	PASS	AA=T;DP=167	GT:GQ:DP	0|1:100:56	0|0:98:33	0|0:100:39
        int start = 201813161 - 1;
        int end = start + 2;
        int count = 0;
        Iterator<Feature> iterator = source.getFeatures("1", start, end);
        while (iterator.hasNext()) {
            VCFVariant feat = (VCFVariant) iterator.next();
            assertEquals("1", feat.getChr());
            assertEquals(start, feat.getStart());
            assertEquals("rs2248941", feat.getID());
            assertEquals(start + 1, feat.getEnd());
            count++;
        }
        assertEquals(1, count);

    }
}