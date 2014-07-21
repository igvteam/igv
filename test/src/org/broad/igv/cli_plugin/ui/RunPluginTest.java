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
