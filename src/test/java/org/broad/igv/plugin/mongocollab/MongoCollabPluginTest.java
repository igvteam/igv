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

package org.broad.igv.plugin.mongocollab;

import com.mongodb.*;
import org.apache.log4j.Logger;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Test;

import java.awt.*;
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

    private static String testFileSpecPath = TestUtils.DATA_DIR + "localMongo.db.spec";
    private static Thread mongoProc;

    private static String MONGO_EXEC_KEY = "MONGO_EXEC_PATH";
    /**
     * See build.xml for default value
     */
    private static String MONGO_EXEC_PATH;

    @BeforeClass
    public static void setUpClass() throws Exception{
        AbstractHeadlessTest.setUpClass();
        MONGO_EXEC_PATH = System.getProperty(MONGO_EXEC_KEY, MONGO_EXEC_PATH);
        log.info("Mongo exec path: " + MONGO_EXEC_PATH);
        Assume.assumeTrue(MONGO_EXEC_PATH != null && MONGO_EXEC_PATH.length() > 0);
        startTestMongo();
        assumeTestDBRunning();
        emptyTestCollection();
    }

    @AfterClass
    public static void tearDownClass() throws Exception{
        AbstractHeadlessTest.tearDownClass();
        try{
            emptyTestCollection();
        }catch(Exception e){
            log.error(e.getMessage(), e);
        }
        try {
            stopTestMongo();
        } catch (Exception e) {
            log.error(e.getMessage(), e);
        }
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
        BasicFeature feat = new BasicFeature("chromo", 50, 100);
        //Set name/desc which look like colors, to make sure color parsing isn't activated
        feat.setColor(Color.magenta);
        feat.setName("0,1,2");
        feat.setDescription("mydescription,is,here");
        DBFeature dbFeat = DBFeature.create(feat);

        assertNull(dbFeat.get_id());

        DBCollection coll = MongoCollabPlugin.getCollection(locator);
        coll.setObjectClass(DBFeature.class);
        String err = MongoCollabPlugin.saveFeature(coll, dbFeat);

        assertNull(err);
        assertNotNull(dbFeat.get_id());

        DBFeature retrievedFeat = (DBFeature) coll.find(new BasicDBObject("_id", dbFeat.get_id())).next();

        TestUtils.assertFeaturesEqual(dbFeat, retrievedFeat);
        assertEquals(dbFeat.getColor(), retrievedFeat.getColor());
        assertEquals(dbFeat.getName(), retrievedFeat.getName());
        assertEquals(dbFeat.getDescription(), retrievedFeat.getDescription());
    }

    @Test
    public void testInsertFeatures() throws Exception{
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

        String dbRelPath = TestUtils.DATA_DIR + "mongodbtmp";
        final File dbDir = new File(dbRelPath);
        String dbAbsPath = dbDir.getAbsolutePath();
        if(!dbDir.exists()){
            dbDir.mkdirs();
            Runtime.getRuntime().addShutdownHook(new Thread(){
                @Override
                public void run() {
                    FileUtils.deleteDir(dbDir);
                }
            });
        }

        final String[] commands = new String[]{mongoExecPath, "--port", "" + port, "--dbpath", dbAbsPath};

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
        //Need to wait for startup
        Thread.sleep(5000);
    }

//    @Test
//    public void testShutdown() throws Exception{
//        assertTrue(canAccessDB(getTestLocator()));
//        assertFalse(stopTestMongo());
//    }


    /**
     * Stop the test Mongo instance
     * @return state of server when done. False = stopped, true = running
     * @throws IOException
     */
    static boolean stopTestMongo() throws Exception{

        MongoCollabPlugin.Locator testLocator = getTestLocator();
        Mongo mongo = MongoCollabPlugin.getMongo(testLocator.host, testLocator.port);
        DB db = mongo.getDB("admin");
        //This actually works (from localhost only), but the java driver throws an error for some reason
        try{
            CommandResult command = db.command( new BasicDBObject( "shutdown", 1 ) );
        }catch(Exception e){

        }
        Thread.sleep(5000);

        if(mongoProc != null && mongoProc.isAlive()){
            mongoProc.interrupt();
            mongoProc.stop();
        }

        return canAccessDB(testLocator);
    }
}
