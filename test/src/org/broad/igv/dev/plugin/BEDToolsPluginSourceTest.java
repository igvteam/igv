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

package org.broad.igv.dev.plugin;

import org.broad.igv.Globals;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.BeforeClass;
import org.junit.Test;
import org.w3c.dom.Element;

import java.io.*;
import java.util.*;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012/05/01
 */
public class BEDToolsPluginSourceTest extends AbstractPluginTest {

    @BeforeClass
    public static void setUpClass() throws Exception {
        pluginPath = "resources/bedtools_plugin.xml";
        AbstractPluginTest.setUpClass();
    }

    @Test
    public void testBedToolsPath() throws Exception {
        String resp = RuntimeUtils.executeShellCommand(new String[]{toolPath}, null, null);
        String line0 = resp.split("\n")[0];
        assertEquals("bedtools: flexible tools for genome arithmetic and DNA sequence analysis.", line0.trim());

    }

    @Test
    public void testGetPath() throws Exception {
        String path = System.getenv("PATH");
        //System.out.println("System path: " + path);
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

    private List<Feature> tstOperationBED3(String cmdArg, int expectedNumFeatures) throws Exception {
        String pathA = TestUtils.DATA_DIR + "bed/test.bed";
        String pathB = TestUtils.DATA_DIR + "bed/test2.bed";
        return tstOperationBED(new String[]{pathA, pathB}, cmdArg, expectedNumFeatures);
    }

    private List<Feature> tstOperationBED6(String cmdArg, int expectedNumFeatures) throws Exception {
        String pathA = TestUtils.DATA_DIR + "bed/H3K4me1_sample_bed6.bed";
        String pathB = TestUtils.DATA_DIR + "bed/H3K4me3_sample_bed6.bed";
        return tstOperationBED(new String[]{pathA, pathB}, cmdArg, expectedNumFeatures);
    }

    private List<Feature> tstOperationBED(String[] paths,
                                          String cmd, int expectedNumFeatures) throws Exception {

        //Find the command element
        Element command = null;
        for (Element curCmd : reader.getCommands(tool)) {
            if (curCmd.getAttribute("cmd").equals(cmd)) {
                command = curCmd;
                break;
            }
        }

        List<Argument> argumentList = reader.getArguments(tool, command);
        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());
        int argnum = 0;
        arguments.put(argumentList.get(argnum), argumentList.get(argnum).getDefaultValue());
        argnum++;

        List<FeatureTrack> trackList = new ArrayList<FeatureTrack>(paths.length);
        TrackLoader loader = new TrackLoader();
        for (String path : paths) {
            TestUtils.createIndex(path);
            FeatureTrack track = (FeatureTrack) loader.load(new ResourceLocator(path), genome).get(0);
            trackList.add(track);
        }

        if (cmd.equals("multiinter")) {
            arguments.put(argumentList.get(argnum++), trackList);
        } else {
            for (FeatureTrack ft : trackList) {
                arguments.put(argumentList.get(argnum++), ft);
            }
        }

        List<String> fullCmd = Arrays.asList(toolPath, cmd);
        PluginFeatureSource combinedFeatureSource = new PluginFeatureSource(fullCmd, arguments,
                reader.getParsingAttributes(tool, command), pluginPath);
        Iterator<Feature> features = combinedFeatureSource.getFeatures("chr1", 0, (int) 1e6);
        List<Feature> featureList = new ArrayList(10);

        while (features.hasNext()) {
            featureList.add(features.next());
        }
        assertEquals("Unexpected number of features in result of " + cmd, expectedNumFeatures, featureList.size());
        return featureList;
    }

    @Test
    public void testOperationsBED3() throws Exception {
        Map<String, Integer> expectedNumFeatures = new HashMap(6);
        expectedNumFeatures.put("intersect", 4);
//        expectedNumFeatures.put("subtract", 2);
//        expectedNumFeatures.put("closest", 6);
//        expectedNumFeatures.put("window", 9);
//        expectedNumFeatures.put("coverage", 3);
//        expectedNumFeatures.put("multiinter", 3);

        for (Map.Entry<String, Integer> entry : expectedNumFeatures.entrySet()) {
            tstOperationBED3(entry.getKey(), entry.getValue());
        }
    }

    @Test
    public void testOperationsBED6() throws Exception {
        Map<String, Integer> expectedNumFeatures = new HashMap(6);
        expectedNumFeatures.put("intersect", 5);
        expectedNumFeatures.put("subtract", 18);
        expectedNumFeatures.put("closest", 19);
        expectedNumFeatures.put("window", 12);
        expectedNumFeatures.put("coverage", 14);
        expectedNumFeatures.put("multiinter", 3);

        for (Map.Entry<String, Integer> entry : expectedNumFeatures.entrySet()) {
            tstOperationBED6(entry.getKey(), entry.getValue());
        }
    }


    @Test
    public void testIntersectBED() throws Exception {
        List<Feature> actFeatures = tstOperationBED3("intersect", 4);
        String expectedPath = TestUtils.DATA_DIR + "bed/isect_res.bed";
        TestUtils.createIndex(expectedPath);
        FeatureSource expFeatureSource = new TribbleFeatureSource(expectedPath, genome);
        Iterator<Feature> expFeatures = expFeatureSource.getFeatures("chr1", 0, (int) 1e6);
        TestUtils.assertFeatureListsEqual(expFeatures, actFeatures.iterator());
    }

    @Test
    public void testClosestBED() throws Exception {
        List<Feature> actFeatures = tstOperationBED3("closest", 6);
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
        String beddir = TestUtils.DATA_DIR + "bed/";
        String[] bedfiles = new String[]{beddir + "test.bed", beddir + "test2.bed", beddir + "isect_res.bed"};

        Iterator<Feature> features = tstOperationBED(bedfiles, "multiinter", 3).iterator();
        int ind = 0;
        int[] expStarts = new int[]{100, 200, 300};
        int[] expEnds = new int[]{101, 201, 301};
        while (features.hasNext()) {
            Feature feat = features.next();
            assertEquals(expStarts[ind], feat.getStart());
            assertEquals(expEnds[ind], feat.getEnd());
            ind++;
        }
    }

    /**
     * Test to make sure we are parsing the 2nd files features
     *
     * @throws Exception
     */
    @Test
    public void testWindow() throws Exception {
        List<Feature> features = tstOperationBED3("window", 9);

        Set<Integer> startsSet = new HashSet<Integer>(Arrays.asList(100, 150, 175));
        Set<Integer> endsSet = new HashSet<Integer>(Arrays.asList(101, 201, 375));
        for (Feature feat : features) {
            assertTrue("Feature start not found in expected set", startsSet.contains(feat.getStart()));
            assertTrue("Feature end not found in expected set", endsSet.contains(feat.getEnd()));
        }
    }

}
