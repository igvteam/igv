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

package org.broad.igv.plugin.mongocollab;

import com.google.common.collect.Lists;
import com.mongodb.DBCollection;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Starts a Mongo server instance, host must have mongo installed
 * See {@link MongoCollabPluginTest}
 * @author jacob
 * @date 2013-Sep-17
 */
public class MongoFeatureSourceTest extends AbstractHeadlessTest{

    private MongoCollabPlugin.Locator locator;
    private MongoFeatureSource source;
    private DBCollection collection;

    @BeforeClass
    public static void setUpClass() throws Exception{
        MongoCollabPluginTest.setUpClass();
    }

    @AfterClass
    public static void tearDownClass() throws Exception{
        MongoCollabPluginTest.tearDownClass();
    }

    @Before
    public void setUp() throws Exception {
        MongoCollabPluginTest.assumeTestDBRunning();
        this.locator = MongoCollabPluginTest.getTestLocator();
        this.collection = MongoCollabPluginTest.emptyTestCollection();
        this.source = new MongoFeatureSource(this.collection, true);
    }

    @After
    public void tearDown() throws Exception{
        MongoCollabPlugin.closeMongo(locator.host, locator.port);
    }

    @Test
    public void testHasIndex() throws Exception{
        assertTrue(this.source.hasIndex());

        this.source = new MongoFeatureSource(this.collection, false);
        assertTrue(this.source.hasIndex());

        this.collection.dropIndexes();
        this.source = new MongoFeatureSource(this.collection, false);
        assertFalse(this.source.hasIndex());
    }

    @Test
    public void testGetFeatures_chr() throws Exception{
        int inserted = MongoCollabPlugin.insertFeaturesFromFile(this.collection, TestUtils.DATA_DIR + "bed/test.bed");

        Iterator<IGVFeature> features = this.source.getFeatures("chr1", 0, Integer.MAX_VALUE);
        List<IGVFeature> chr1List = Lists.newArrayList(features);

        features = this.source.getFeatures("chr2", 0, Integer.MAX_VALUE);
        List<IGVFeature> chr2List = Lists.newArrayList(features);

        assertEquals(inserted, chr1List.size() + chr2List.size());

        TestUtils.assertFeatureIteratorSorted(chr1List.iterator());

        TestUtils.assertFeatureIteratorSorted(chr2List.iterator());
    }

    @Test
    public void testGetFeatures_start_01() throws Exception{
        int inserted = MongoCollabPlugin.insertFeaturesFromFile(this.collection, TestUtils.DATA_DIR + "bed/test.bed");

        Iterator<IGVFeature> features = this.source.getFeatures("chr1", 250, 100005);
        List<IGVFeature> list = Lists.newArrayList(features);

        assertEquals(2, list.size());
        TestUtils.assertFeatureIteratorSorted(list.iterator());
    }

    @Test
    public void testGetFeatures_start_02() throws Exception{
        int inserted = MongoCollabPlugin.insertFeaturesFromFile(this.collection, TestUtils.DATA_DIR + "bed/test.bed");

        Iterator<IGVFeature> features = this.source.getFeatures("chr1", 100005, 100008);
        List<IGVFeature> list_00 = Lists.newArrayList(features);

        assertEquals(1, list_00.size());
        assertEquals(100000, list_00.get(0).getStart());

        features = this.source.getFeatures("chr1", 100005, 200008);
        List<IGVFeature> list_01 = Lists.newArrayList(features);

        assertEquals(2, list_01.size());
        assertEquals(100000, list_01.get(0).getStart());
    }
}
