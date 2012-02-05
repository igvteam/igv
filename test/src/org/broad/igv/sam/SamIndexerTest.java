/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author jrobinso
 */
public class SamIndexerTest {

    public SamIndexerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of createSamIndex method, of class SamIndexer.
     */
    //@Test
    public void createSamIndex() throws Exception {
        /*
        System.out.println("createSamIndex");
        File samFile = new File("/Users/jrobinso/IGV/Sam/test.sam");
        File idxFile = new File("/Users/jrobinso/IGV/Sam/test.sai");
        int tileWidth = 100;
        AlignmentIndexer instance = AlignmentIndexer.getInstance(samFile, null, null);
        FeatureIndex result = instance.createSamIndex(idxFile, tileWidth);
        // TODO review the generated test code and remove the default call to fail.

        result.store(idxFile);


        SamQueryTextReader reader = new SamQueryTextReader(samFile, idxFile);

        CloseableIterator<Alignment> iter = reader.query("chrM", 0, 200, false);

        while (iter.hasNext()) {
            Alignment rec = iter.next();
            System.out.println(rec.getReadName());
        }

        iter.close();
        */

    }


    @Test
    public void test2() throws Exception {
        /*
        System.out.println("createSamIndex");
        File samFile = new File("/Users/jrobinso/IGV/Sam/303KY.8.paired.sam");
        File idxFile = new File("/Users/jrobinso/IGV/Sam/303KY.8.paired.sam.sai");
        int tileWidth = 10000;
        AlignmentIndexer instance = AlignmentIndexer.getInstance(samFile, null, null);
        FeatureIndex result = instance.createSamIndex(idxFile, tileWidth);
        result.store(idxFile);


        SamQueryTextReader reader = new SamQueryTextReader(samFile, idxFile);

        CloseableIterator<Alignment> iter = reader.query("chr1", 100, 200, false);

        while (iter.hasNext()) {
            Alignment rec = iter.next();
            System.out.println(rec.getReadName());
        }

        iter.close();
        */

    }


}
