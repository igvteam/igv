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

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.BeforeClass;
import org.junit.Test;
import org.w3c.dom.Element;

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
        Element command = null;
        for (Element curCmd : reader.getCommands(tool)) {
            if (curCmd.getAttribute("name").equals("Filter by size")) {
                command = curCmd;
                break;
            }
        }

        Map<String, String> parsingAttrs = reader.getParsingAttributes(tool, command);

        String testFile = TestUtils.DATA_DIR + "bed/Unigene.sample.bed";
        FeatureTrack track = (FeatureTrack) (new TrackLoader()).load(new ResourceLocator(testFile), genome).get(0);

        List<Argument> argumentList = reader.getArguments(tool, command);
        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());
        int argnum = 0;

        Object[] argVals = {"" + min, "" + max, track};
        for (Object arg : argVals) {
            arguments.put(argumentList.get(argnum), arg);
            argnum++;
        }

        List<String> fullCmd = Arrays.asList(toolPath);
        PluginFeatureSource source = new PluginFeatureSource(fullCmd, arguments, parsingAttrs, pluginPath);
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
