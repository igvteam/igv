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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.*;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012/05/01
 */
@Ignore
public class CombinedFeatureSourceTest extends AbstractHeadlessTest {

    static boolean haveBedTools;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        haveBedTools = CombinedFeatureSource.checkBEDToolsPathValid();
        Assume.assumeTrue(haveBedTools);
    }

    @Test
    public void testBedToolsPath() throws Exception {
        String cmd = FileUtils.findExecutableOnPath(Globals.BEDtoolsPath);
        String resp = RuntimeUtils.executeShellCommand(new String[]{cmd}, null, null);
        String line0 = resp.split("\n")[0];
        assertEquals("bedtools: flexible tools for genome arithmetic and DNA sequence analysis.", line0.trim());

    }

    @Test
    public void testGetPath() throws Exception {
        String path = System.getenv("PATH");
        System.out.println("System path: " + path);
        assertNotNull(path);
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
            System.out.println(line);
            errlines++;
        }

        assertEquals(2, readlines);
        assertEquals(0, errlines);

    }

    private List<Feature> tstOperationBED3(CombinedFeatureSource.Operation operation, int expectedNumFeatures) throws Exception {
        String pathA = TestUtils.DATA_DIR + "bed/test.bed";
        String pathB = TestUtils.DATA_DIR + "bed/test2.bed";
        return tstOperationBED(pathA, pathB, operation, expectedNumFeatures);
    }

    private List<Feature> tstOperationBED6(CombinedFeatureSource.Operation operation, int expectedNumFeatures) throws Exception {
        String pathA = TestUtils.DATA_DIR + "bed/H3K4me1_sample_bed6.bed";
        String pathB = TestUtils.DATA_DIR + "bed/H3K4me3_sample_bed6.bed";
        return tstOperationBED(pathA, pathB, operation, expectedNumFeatures);
    }

    private List<Feature> tstOperationBED(String pathA, String pathB,
                                          CombinedFeatureSource.Operation operation, int expectedNumFeatures) throws Exception {
        IgvTools igvTools = new IgvTools();
        igvTools.doIndex(pathA, null, 1, 16000);
        igvTools.doIndex(pathB, null, 1, 16000);
        FeatureSource sourceA = TribbleFeatureSource.getFeatureSource(new ResourceLocator(pathA), genome);
        FeatureSource sourceB = TribbleFeatureSource.getFeatureSource(new ResourceLocator(pathB), genome);


        CombinedFeatureSource combinedFeatureSource = new CombinedFeatureSource(new FeatureSource[]{sourceA, sourceB}, operation);
        Iterator<Feature> features = combinedFeatureSource.getFeatures("chr1", 0, (int) 1e6);
        List<Feature> featureList = new ArrayList(10);

        while (features.hasNext()) {
            featureList.add(features.next());
        }
        assertEquals("Unexpected number of features in result of " + operation, expectedNumFeatures, featureList.size());
        return featureList;
    }

    @Test
    public void testOperationsBED3() throws Exception {
        Map<CombinedFeatureSource.Operation, Integer> expectedNumFeatures = new HashMap(6);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.INTERSECT, 4);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.SUBTRACT, 2);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.CLOSEST, 6);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.WINDOW, 9);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.COVERAGE, 3);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.MULTIINTER, 3);

        for (Map.Entry<CombinedFeatureSource.Operation, Integer> entry : expectedNumFeatures.entrySet()) {
            tstOperationBED3(entry.getKey(), entry.getValue());
        }
    }

    @Test
    public void testOperationsBED6() throws Exception {
        Map<CombinedFeatureSource.Operation, Integer> expectedNumFeatures = new HashMap(6);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.INTERSECT, 5);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.SUBTRACT, 18);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.CLOSEST, 19);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.WINDOW, 12);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.COVERAGE, 14);
        expectedNumFeatures.put(CombinedFeatureSource.Operation.MULTIINTER, 3);

        for (Map.Entry<CombinedFeatureSource.Operation, Integer> entry : expectedNumFeatures.entrySet()) {
            tstOperationBED6(entry.getKey(), entry.getValue());
        }
    }


    @Test
    public void testIntersectBED() throws Exception {
        List<Feature> actFeatures = tstOperationBED3(CombinedFeatureSource.Operation.INTERSECT, 4);
        String expectedPath = TestUtils.DATA_DIR + "bed/isect_res.bed";
        FeatureSource expFeatureSource = TribbleFeatureSource.getFeatureSource(new ResourceLocator(expectedPath), genome);
        Iterator<Feature> expFeatures = expFeatureSource.getFeatures("chr1", 0, (int) 1e6);
        TestUtils.assertFeatureListsEqual(expFeatures, actFeatures.iterator());

    }

    @Test
    public void testClosestBED() throws Exception {
        List<Feature> actFeatures = tstOperationBED3(CombinedFeatureSource.Operation.CLOSEST, 6);
        String expectedPath = TestUtils.DATA_DIR + "bed/closest_res.bed";
        BufferedReader reader = new BufferedReader(new FileReader(expectedPath));
        int ind = 0;
        String line = null;
        while ((line = reader.readLine()) != null) {

            String[] tokens = Globals.singleTabMultiSpacePattern.split(line);

            String expChr = tokens[3];
            int expStart = Integer.parseInt(tokens[4]);
            int expEnd = Integer.parseInt(tokens[5]);

            if (expChr.equalsIgnoreCase(".")) {
                continue;
            }

            Feature actFeature = actFeatures.get(ind);
            assertEquals(expChr, actFeature.getChr());
            assertEquals(expStart, actFeature.getStart());
            assertEquals(expEnd, actFeature.getEnd());

            ind++;

        }

    }

    @Test
    public void testMultiinterBED() throws Exception {
        String[] bedfiles = new String[]{"test.bed", "test2.bed", "isect_res.bed"};
        List<FeatureSource> sources = new ArrayList<FeatureSource>(bedfiles.length);
        for (String bedfile : bedfiles) {
            sources.add(TribbleFeatureSource.getFeatureSource(new ResourceLocator(TestUtils.DATA_DIR + "bed/" + bedfile), genome));
        }

        CombinedFeatureSource combinedFeatureSource = new CombinedFeatureSource(sources.toArray(new FeatureSource[0]),
                CombinedFeatureSource.Operation.MULTIINTER);
        Iterator<Feature> features = combinedFeatureSource.getFeatures("chr1", 0, (int) 1e6);
        int ind = 0;
        int[] expStarts = new int[]{100, 200, 300};
        int[] expEnds = new int[]{101, 201, 301};
        while (features.hasNext()) {
            Feature feat = features.next();
            assertEquals(expStarts[ind], feat.getStart());
            assertEquals(expEnds[ind], feat.getEnd());
            ind++;
        }

        assertEquals(3, ind);

    }

    /**
     * Test to make sure we are parsing the 2nd files features
     *
     * @throws Exception
     */
    @Test
    public void testWindow() throws Exception {
        List<Feature> features = tstOperationBED3(CombinedFeatureSource.Operation.WINDOW, 9);

        Set<Integer> startsSet = new HashSet<Integer>(Arrays.asList(100, 150, 175));
        Set<Integer> endsSet = new HashSet<Integer>(Arrays.asList(101, 201, 375));
        for (Feature feat : features) {
            assertTrue("Feature start not found in expected set", startsSet.contains(feat.getStart()));
            assertTrue("Feature end not found in expected set", endsSet.contains(feat.getEnd()));
        }
    }

}
