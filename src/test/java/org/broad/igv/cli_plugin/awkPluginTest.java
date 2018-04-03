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

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.*;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Aug-20
 */
public class awkPluginTest extends AbstractPluginTest {

    @BeforeClass
    public static void setUpClass() throws Exception {
        pluginPath = "resources/awk_plugin.xml";
        AbstractPluginTest.setUpClass();
    }

    @Test
    public void testFilterBySize() throws Exception {
        int[] mins = {0, 1, 50, 200000000 - 1};
        int[] maxs = {0, 500000, 5000, 200000000 - 1};
        int[] expNums = {0, 71, 57, 1};
        for (int ind = 0; ind < mins.length; ind++) {
            tstFilterSize(mins[ind], maxs[ind], expNums[ind]);
        }

    }

    public List<Feature> tstFilterSize(int min, int max, int expectedNum) throws Exception {

        //Find the command element
        PluginSpecReader.Command  command = null;
        for (PluginSpecReader.Command  curCmd : tool.commandList) {
            if (curCmd.name.equals("Filter by size")) {
                command = curCmd;
                break;
            }
        }

        List<PluginSpecReader.Output> outputAttrs = command.outputList;
        PluginSpecReader.Output outputAttr = outputAttrs.get(0);

        String testFile = TestUtils.DATA_DIR + "bed/Unigene.sample.bed";
        FeatureTrack track = (FeatureTrack) (new TrackLoader()).load(new ResourceLocator(testFile), genome).get(0);

        List<Argument> argumentList = command.argumentList;
        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());
        int argnum = 0;

        Object[] argVals = {"" + min, "" + max, track};
        for (Object arg : argVals) {
            arguments.put(argumentList.get(argnum), arg);
            argnum++;
        }

        List<String> fullCmd = Arrays.asList(toolPath);
        PluginFeatureSource source = new PluginFeatureSource(fullCmd, arguments, outputAttr, pluginPath);
        List<Feature> feats = new ArrayList<Feature>();
        Iterator<Feature> featIter = source.getFeatures("chr2", 0, Integer.MAX_VALUE);
        Feature feat;
        while (featIter.hasNext()) {
            feat = featIter.next();
            int length = feat.getEnd() - feat.getStart();
            assertTrue(length >= min);
            assertTrue(length <= max);
            feats.add(feat);
        }
        assertEquals(expectedNum, feats.size());
        return feats;
    }


}
