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
import org.apache.log4j.Logger;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.*;

import java.io.*;

import static org.junit.Assert.*;

/**
 * Tests of our integration with MongoDB for collaboration/annotation
 * Starts MongoDB on a separate process, host machine must have Mongo and
 * test runner must be pointed to it by {@link #MONGO_EXEC_KEY}
 * @author jacob
 * @date 2013-Sep-10
 */
public class MongoCollabPluginTest extends AbstractHeadlessTest {

    private static Logger log = Logger.getLogger(MongoCollabPluginTest.class);

    private static String testFileSpecPath = TestUtils.DATA_DIR + "testMongo.db.spec";
    private static Thread mongoProc;

    private static String MONGO_EXEC_KEY = "MONGO_EXEC_PATH";
    /**
     * See build.xml for default value
     */
    private static String MONGO_EXEC_PATH;

    @BeforeClass
    public static void setUpClass() throws Exception{
        MONGO_EXEC_PATH = System.getProperty(MONGO_EXEC_KEY, MONGO_EXEC_PATH);
        log.info("Mongo exec path: " + MONGO_EXEC_PATH);
        Assume.assumeTrue(MONGO_EXEC_PATH != null && MONGO_EXEC_PATH.length() > 0);
        startTestMongo();
        assumeTestDBRunning();
    }

    @AfterClass
    public static void tearDownClass() throws Exception{
        stopTestMongo();
        TestUtils.clearOutputDir();
    }

    @Before
    public void setUp() throws Exception{
        assumeTestDBRunning();
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
        DBFeature dbFeat = DBFeature.create(feat);

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

    public static MongoCollabPlugin.Locator getTestLocator() throws IOException{
        try {
            return new MongoCollabPlugin.Locator(testFileSpecPath);
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

    static void startTestMongo() throws Exception{
        String mongoExecPath = MONGO_EXEC_PATH;
        int port = getTestLocator().port;
        String dbPath = TestUtils.TMP_OUTPUT_DIR;
        File dbDir = new File(dbPath);
        if(!dbDir.exists()){
            dbDir.mkdirs();
            dbDir.deleteOnExit();
        }

        final String[] commands = new String[]{mongoExecPath, "--port", "" + port, "--dbpath", dbPath};

        Thread runnable = new Thread(){

            private Process proc;

            @Override
            public void run() {
                try {
                    proc = RuntimeUtils.startExternalProcess(commands, null, null);
                    InputStream is = proc.getInputStream();
                    BufferedReader br = new BufferedReader(new InputStreamReader(is));
                    String line = null;
                    while((line = br.readLine()) != null){
                        log.debug(line);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            @Override
            public void interrupt() {
                if(this.proc != null){
                    this.proc.destroy();
                }
                super.interrupt();
            }
        };
        runnable.start();
        mongoProc = runnable;
    }


    static void stopTestMongo() {
        if(mongoProc != null){
            mongoProc.interrupt();
            mongoProc.stop();
        }
    }
}
