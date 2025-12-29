package org.broad.igv.tools.motiffinder;

import org.broad.igv.track.Track;
import org.broad.igv.ui.AbstractHeadedTest;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Oct-11
 */
public class MotifFinderPluginHeadlessTest extends AbstractHeadedTest{

    @Test
    public void addMultiTracks() throws Exception {
        String[] patterns = new String[]{"TATTAAATTA", "CGCGCGCGCCGNATNG", "N.[A,C]CNCGCTC"};
        String[] posNames = new String[]{"tata", "cg", "regex"};
        String[] negNames = new String[]{"tata neg", "cg neg", "regex neg"};
        List<Track> tracks = MotifFinderPlugin.generateTracksForPatterns(patterns, posNames, negNames);
        assertEquals(2*patterns.length, tracks.size());
    }
}
