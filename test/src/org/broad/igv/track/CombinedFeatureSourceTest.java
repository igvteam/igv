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

package org.broad.igv.track;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;

import static org.junit.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/05/01
 */
public class CombinedFeatureSourceTest {

    Genome genome;

    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();
        genome = TestUtils.loadGenome();

    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testBedToolsPath() throws Exception {
        String cmd = CombinedFeatureSource.BEDtoolsPath;
        String resp = RuntimeUtils.executeShellCommand(cmd, null, null);
        String line0 = resp.split("\n")[0];
        assertEquals("bedtools: flexible tools for genome arithmetic and DNA sequence analysis.", line0.trim());

    }

    //Test our ability to write to stdin and read from stdout
    //XXX This test won't work on windows
    //TODO Move this to RuntimeUtilsTest (which doesn't exist yet)
    @Test
    public void testWriteStdIn() throws Exception {
        if (Globals.IS_WINDOWS) {
            return;
        }
        Runtime run = Runtime.getRuntime();
        Process pr = null;

        String msg = "Never have I ever \n thought twas ever thus \n a good movie \n thus I'm glad it's over";
        try {
            pr = run.exec("/usr/bin/grep -i thus", null, null);
        } catch (IOException e) {
            e.printStackTrace();
        }

        BufferedReader in = new BufferedReader(new InputStreamReader(pr.getInputStream()));
        BufferedReader err = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(pr.getOutputStream())), true);
        out.println(msg);

        String line;

        out.flush();
        out.close();
        pr.waitFor();

        int readlines = 0;
        while ((line = in.readLine()) != null) {
            //System.out.println(line);
            readlines++;
        }

        in.close();


        //System.out.println("errors:");

        int errlines = 0;
        while ((line = err.readLine()) != null) {
            //System.out.println(line);
            errlines++;
        }

        assertEquals(2, readlines);
        assertEquals(0, errlines);

    }

    @Ignore
    @Test
    public void testIntersect() throws Exception {

        String pathA = TestUtils.DATA_DIR + "bed/test.bed";
        String pathB = TestUtils.DATA_DIR + "bed/test2.bed";
        FeatureSource sourceA = new TribbleFeatureSource(pathA, genome);
        FeatureSource sourceB = new TribbleFeatureSource(pathB, genome);

        CombinedFeatureSource.Operation operation = CombinedFeatureSource.Operation.INTERSECT;

        CombinedFeatureSource combinedFeatureSource = new CombinedFeatureSource(sourceA, sourceB, operation);


    }
}
