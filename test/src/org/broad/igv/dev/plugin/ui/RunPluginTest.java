package org.broad.igv.dev.plugin.ui;

import org.broad.igv.dev.plugin.AbstractPluginTest;
import org.broad.igv.dev.plugin.PluginSpecReader;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.tribble.Feature;
import org.junit.Ignore;
import org.junit.Test;
import org.w3c.dom.Element;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**.
 * User: jacob
 * Date: Nov 26 2012
 */
@Ignore
public class RunPluginTest extends AbstractHeadedTest {

    @Test
    public void testCatPlugin() throws Exception {
        PluginSpecReader catReader = AbstractPluginTest.getCatReader();
        Element tool = catReader.getTools().get(0);
        Element command = catReader.getCommands(tool).get(0);
        RunPlugin rp = new RunPlugin(IGV.getMainFrame(), catReader, tool, command);

        int numTracksBefore = IGV.getInstance().getVisibleTrackCount();
        rp.okButtonActionPerformed(null);

        int numTracksAfter = IGV.getInstance().getVisibleTrackCount();
        assertEquals(numTracksBefore + 1, numTracksAfter);

        FeatureTrack newTrack = (FeatureTrack) IGV.getInstance().getAllTracks().get(numTracksAfter-1);
        List<Feature> features = newTrack.getFeatures("chr5", 1, 10000);
        assertNotNull(features);
        assertTrue(features.size() > 0);
    }
}
