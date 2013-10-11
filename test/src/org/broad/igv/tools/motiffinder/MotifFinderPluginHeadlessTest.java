/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

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
