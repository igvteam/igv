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

package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2013-Feb-21
 */
public class FeatureSearcherTest extends AbstractHeadlessTest {

    @Test
    public void testSearchInWindow() throws Exception{
        tstSearchInWindow(true);
        tstSearchInWindow(false);
    }

    @Test
    public void testSearchOutsideWindow() throws Exception{
        tstSearchOutsideWindow(true);
        tstSearchOutsideWindow(false);
    }

    @Test
    public void testSearchOutsideWindowBackwards() throws Exception{
        //tstSearchOutsideWindowBackwards(true);
        tstSearchOutsideWindowBackwards(false);
    }

    @Test
    public void testSearchNextChromo() throws Exception{
        tstSearchDifferentChromo(true);
        tstSearchDifferentChromo(false);
    }

    /**
     * Runs the search, either on this thread or a different one.
     * In either case blocks until the searching is done
     * @param searcher
     * @param sepThread
     * @throws Exception
     */
    private void runSearch(FeatureSearcher searcher, boolean sepThread) throws Exception{

        assertFalse(searcher.isDone());
        assertNull(searcher.getResult());

        if(sepThread){
            LongRunningTask.getThreadExecutor().execute(searcher);
            //Tacky, but it should take <1000 mSec to start the search
            Thread.sleep(1000);
            while(!searcher.isDone()){
                Thread.sleep(100);
            }
        }else{
            searcher.run();
        }
    }

    public void tstSearchInWindow(boolean sepThread) throws Exception{

        FeatureSearcher searcher = new FeatureSearcher(getTestBedSource(), genome, "chr1", 0);

        runSearch(searcher, sepThread);

        Feature feat = searcher.getResult().next();

        assertEquals("chr1", feat.getChr());
        assertEquals(100, feat.getStart());
        assertEquals(101, feat.getEnd());
    }

    public void tstSearchOutsideWindow(boolean sepThread) throws Exception{

        FeatureSearcher searcher = new FeatureSearcher(getTestBedSource(), genome, "chr1", 500);
        searcher.setSearchIncrement(10000);

        runSearch(searcher, sepThread);

        Feature feat = searcher.getResult().next();

        assertEquals("chr1", feat.getChr());
        assertEquals(100000, feat.getStart());
        assertEquals(100010, feat.getEnd());
    }

    public void tstSearchOutsideWindowBackwards(boolean sepThread) throws Exception{

        FeatureSearcher searcher = new FeatureSearcher(getTestBedSource(), genome, "chr1", 500000);
        searcher.setSearchIncrement(-10000);

        runSearch(searcher, sepThread);

        Feature feat = searcher.getResult().next();

        assertEquals("chr1", feat.getChr());
        assertEquals(100000, feat.getStart());
        assertEquals(100010, feat.getEnd());
    }

    public void tstSearchDifferentChromo(boolean sepThread) throws Exception{

        FeatureSearcher searcher = new FeatureSearcher(getTestBedSource(), genome, "chr1", 500000);
        searcher.setSearchIncrement(10000);

        runSearch(searcher, sepThread);

        Feature feat = searcher.getResult().next();

        assertEquals("chr2", feat.getChr());
        assertEquals(1, feat.getStart());
        assertEquals(10, feat.getEnd());
    }

    public FeatureSource<? extends Feature> getTestBedSource() throws IOException{
        String path = TestUtils.DATA_DIR + "bed/test.bed";
        TestUtils.createIndex(path);
        return new TribbleFeatureSource(path, genome);
    }

    /**
     * Test that stopping actually stops the search.
     * We start searching a window which will never complete
     * @throws Exception
     */
    @Test
    public void testCancel() throws Exception{
        FeatureSearcher searcher = new FeatureSearcher(getTestBedSource(), genome, "chr3", 0);
        searcher.setSearchIncrement(100);

        // Be careful monkeying around here.
        // The test will exit and shut down the threadExecutor if there is nothing blocking on this thread
        // This can make it seem like the search aborted prematurely
        LongRunningTask.getThreadExecutor().execute(searcher);
        Thread.sleep(500);
        assertFalse(searcher.isDone());
        searcher.cancel();

        //Stopping is not instantaneous
        int counter = 0;
        while(!searcher.isDone()){
            Thread.sleep(100);
            counter++;
        }
        assertNull(searcher.getResult());
        assertTrue(counter < 5);
    }
}
