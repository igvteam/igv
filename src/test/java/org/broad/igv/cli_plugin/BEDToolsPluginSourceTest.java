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

package org.broad.igv.cli_plugin;

import org.broad.igv.Globals;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
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

        PluginSpecReader.Command command = findCommandElementByName(tool, cmd);
        List<Argument> argumentList = command.argumentList;
        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());
        int argnum = 0;
        int numArgsSet = 0;
        if(cmd.equalsIgnoreCase("intersect")){
            numArgsSet = 2;
        }else if(cmd.equalsIgnoreCase("coverage")){
            numArgsSet = 1;
        }
        for(; argnum < numArgsSet; argnum++){

            Argument argument = argumentList.get(argnum);
            Object value = argumentList.get(argnum).getDefaultValue();

            if(argument.getType() == Argument.InputType.BOOL){
                value = Boolean.parseBoolean((String) value);
            }

            arguments.put(argument, value);
        }


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

        //Optional arguments
        arguments.put(argumentList.get(argnum), argumentList.get(argnum).getDefaultValue());

        List<String> commands = Arrays.asList(toolPath, cmd);
        PluginFeatureSource combinedFeatureSource = new PluginFeatureSource(commands, arguments, command.outputList.get(0), pluginPath);
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
        expectedNumFeatures.put("subtract", 2);
        expectedNumFeatures.put("closest", 6);
        expectedNumFeatures.put("window", 9);
        expectedNumFeatures.put("coverage", 3);
        expectedNumFeatures.put("multiinter", 3);

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
        FeatureSource expFeatureSource = TribbleFeatureSource.getFeatureSource(new ResourceLocator(expectedPath), genome);
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
