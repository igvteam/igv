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

package org.broad.igv.cli_plugin.ui;

import org.broad.igv.cli_plugin.AbstractPluginTest;
import org.broad.igv.cli_plugin.PluginSpecReader;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.*;

/**
 * Test our RunPlugin dialog
 * User: jacob
 * Date: Nov 26 2012
 */
public class RunPluginTest extends AbstractHeadedTest {

    @Test
    public void testCatPlugin() throws Exception {
        PluginSpecReader catReader = AbstractPluginTest.getCatReader();
        AbstractPluginTest.initTool(catReader);

        PluginSpecReader.Tool tool = catReader.getTools().get(0);
        PluginSpecReader.Command command = tool.commandList.get(0);
        RunPlugin rp = new RunPlugin(IGV.getMainFrame(), catReader, tool, command);

        int numTracksBefore = IGV.getInstance().getAllTracks().size();
        rp.okButtonActionPerformed(null);

        int numTracksAfter = IGV.getInstance().getAllTracks().size();
        assertEquals(numTracksBefore + 1, numTracksAfter);

        FeatureTrack newTrack = (FeatureTrack) IGV.getInstance().getAllTracks().get(numTracksAfter - 1);
        List<Feature> features = newTrack.getFeatures("chr5", 1, 10000);
        assertNotNull(features);
        assertTrue(features.size() > 0);
    }
}
