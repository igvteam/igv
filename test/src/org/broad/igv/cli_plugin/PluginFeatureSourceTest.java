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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.Assume;
import org.junit.Test;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Aug-10
 */
public class PluginFeatureSourceTest extends AbstractHeadlessTest {

    @Test
    public void testGetFeaturesCat() throws Exception {

        Assume.assumeTrue(!Globals.IS_WINDOWS);

        PluginSpecReader reader = AbstractPluginTest.getCatReader();

        PluginSpecReader.Tool tool = reader.getTools().get(0);
        PluginSpecReader.Command command = tool.commandList.get(0);
        List<Argument> argumentList = command.argumentList;


        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());

        int argnum = 0;
        arguments.put(argumentList.get(argnum++), "");

        TrackLoader loader = new TrackLoader();
        String[] paths = new String[]{TestUtils.DATA_DIR + "bed/test.bed", TestUtils.DATA_DIR + "bed/testAlternateColor.bed"};
        for (String path : paths) {
            TestUtils.createIndex(path);
            FeatureTrack track = (FeatureTrack) loader.load(new ResourceLocator(path), genome).get(0);
            arguments.put(argumentList.get(argnum++), track);
        }


        List<String> cmd = Arrays.asList(reader.getToolPath(tool), command.cmd);
        PluginFeatureSource pluginSource = new PluginFeatureSource(cmd, arguments, command.outputList.get(0), reader.getSpecPath());

        Iterator<Feature> featuresExp = ((FeatureTrack) arguments.get(argumentList.get(1))).getFeatures("chr2", 1, 30).iterator();
        Iterator<Feature> featuresAct = pluginSource.getFeatures("chr2", 1, 30);

        TestUtils.assertFeatureListsEqual(featuresExp, featuresAct);

        int start = 178707289 - 1;
        int end = 179714478;
        featuresExp = ((FeatureTrack) arguments.get(argumentList.get(2))).getFeatures("chr2", start, end).iterator();
        featuresAct = pluginSource.getFeatures("chr2", start, end);

        TestUtils.assertFeatureListsEqual(featuresExp, featuresAct);

    }
}
