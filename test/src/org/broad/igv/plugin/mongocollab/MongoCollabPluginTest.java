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

import com.mongodb.DBCollection;
import com.mongodb.Mongo;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;

import static org.junit.Assert.*;

/**
 * Tests of our integration with MongoDB for collaboration/annotation
 * @author jacob
 * @date 2013-Sep-10
 */
@Ignore("DB must be started manually")
public class MongoCollabPluginTest extends AbstractHeadlessTest {

    private static String testFileSpecPath = TestUtils.DATA_DIR + "testMongo.txt";

    @Before
    public void setUp() throws Exception{
        //assumeTestDBRunning();
    }

    @Test
    public void testReadSpec() throws Exception{
        MongoCollabPlugin.Locator locator = getTestLocator();
        assertEquals("localhost", locator.host);
        assertEquals(27017, locator.port);
    }

    @Test
    public void testGetCollection() throws Exception{
        MongoCollabPlugin.Locator locator = assumeTestDBRunning();
        DBCollection collection = MongoCollabPlugin.getCollection(locator);
        assertNotNull(collection);
    }

    @Test
    public void testInsertFeature() throws Exception{
        MongoCollabPlugin.Locator locator = assumeTestDBRunning();
        Feature feat = new BasicFeature("chromo", 50, 100);
        MongoCollabPlugin.FeatDBObject dbFeat = MongoCollabPlugin.FeatDBObject.create(feat);

        assertNull(dbFeat.get_id());

        String err = MongoCollabPlugin.saveFeature(MongoCollabPlugin.getCollection(locator), dbFeat);

        assertNull(err);
        assertNotNull(dbFeat.get_id());

        TestUtils.assertFeaturesEqual(feat, dbFeat);
    }

    @Test
    public void testInsertFeatures() throws Exception{
        assumeTestDBRunning();
        DBCollection collection = emptyTestCollection();
        int inserted = MongoCollabPlugin.insertFeaturesFromFile(collection, TestUtils.DATA_DIR + "bed/test.bed");

        assertEquals(inserted, collection.count());
        assertTrue(inserted > 0);
    }

    /**
     * Assumes that test DB is running, return Locator (for convenience)
     * if it is.
     * @return
     * @throws Exception
     */
    public static MongoCollabPlugin.Locator assumeTestDBRunning() throws Exception{
        MongoCollabPlugin.Locator locator = getTestLocator();
        Assume.assumeTrue(canAccessDB(locator));
        return locator;
    }

    public static DBCollection emptyTestCollection() throws Exception{
        DBCollection collection = MongoCollabPlugin.getCollection(getTestLocator());
        collection.dropIndexes();
        collection.drop();
        return collection;
    }

    public static MongoCollabPlugin.Locator getTestLocator(){
        try {
            return new MongoCollabPlugin.Locator(new File(testFileSpecPath));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static boolean canAccessDB(MongoCollabPlugin.Locator locator){
        try {
            Mongo mongo = MongoCollabPlugin.getMongo(locator.host, locator.port);
            return mongo.getDatabaseNames() != null;
        } catch (Exception e) {
            return false;
        }
    }
}
