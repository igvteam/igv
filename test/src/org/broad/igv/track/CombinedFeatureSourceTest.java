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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/05/01
 */
@Ignore
public class CombinedFeatureSourceTest extends AbstractHeadlessTest {


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

    public List<Feature> tstOperationBED(CombinedFeatureSource.Operation operation, int expectedNumFeatures) throws Exception {
        String pathA = TestUtils.DATA_DIR + "bed/test.bed";
        String pathB = TestUtils.DATA_DIR + "bed/test2.bed";
        IgvTools igvTools = new IgvTools();
        igvTools.doIndex(pathA, 1, 16000);
        igvTools.doIndex(pathB, 1, 16000);
        FeatureSource sourceA = new TribbleFeatureSource(pathA, genome);
        FeatureSource sourceB = new TribbleFeatureSource(pathB, genome);


        CombinedFeatureSource combinedFeatureSource = new CombinedFeatureSource(sourceA, sourceB, operation);
        Iterator<Feature> features = combinedFeatureSource.getFeatures("chr1", 0, (int) 1e6);
        List<Feature> featureList = new ArrayList(10);

        while (features.hasNext()) {
            featureList.add(features.next());
        }
        assertEquals("Unexpected number of features in result of " + operation, expectedNumFeatures, featureList.size());
        return featureList;
    }

    @Test
    public void testOperationsBED() throws Exception {
        Map<CombinedFeatureSource.Operation, Integer> expectedNumFeatures = new HashMap(4);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.INTERSECT, 4);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.SUBTRACT, 2);
        //expectedNumFeatures.put(CombinedFeatureSource.Operation.CLOSEST, 6);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.WINDOW, 9);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.COVERAGE, 3);

        for (Map.Entry<CombinedFeatureSource.Operation, Integer> entry : expectedNumFeatures.entrySet()) {
            tstOperationBED(entry.getKey(), entry.getValue());
        }
    }

    /**
     * Test to make sure we are parsing the 2nd files features
     *
     * @throws Exception
     */
    @Test
    public void testWindow() throws Exception {
        List<Feature> features = tstOperationBED(CombinedFeatureSource.Operation.WINDOW, 9);

        Set<Integer> startsSet = new HashSet<Integer>(Arrays.asList(100, 150, 175));
        Set<Integer> endsSet = new HashSet<Integer>(Arrays.asList(101, 201, 375));
        for (Feature feat : features) {
            assertTrue("Feature start not found in expected set", startsSet.contains(feat.getStart()));
            assertTrue("Feature end not found in expected set", endsSet.contains(feat.getEnd()));
        }
    }

}
