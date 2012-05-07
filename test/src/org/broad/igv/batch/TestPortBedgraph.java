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

package org.broad.igv.batch;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import static org.junit.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 11/30/11
 */
public class TestPortBedgraph {

    private BufferedReader in;

    public TestPortBedgraph() {
    }


    public int importBedGraph(String filename) throws Exception {
        in = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
        String line = "";
        int count = 0;
        while ((line = in.readLine()) != null) {
            count++;
        }
        return count;
    }

    @Test
    public void test1409() throws Exception {
        String testfile = TestUtils.DATA_DIR + "wig/jira_1409.bedgraph";
        int count = importBedGraph(testfile);
        assertTrue(count > 0);
    }


}
