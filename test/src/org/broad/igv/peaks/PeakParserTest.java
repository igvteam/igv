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

package org.broad.igv.peaks;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 9/28/11
 */
public class PeakParserTest {
    @Test
    public void testLoadPeaksBinary() throws Exception {

        String f = TestUtils.DATA_DIR + "ichip/stat6.peak.bin";
        PeakParser parser = new PeakParser(f);
        List<Peak> peaks = parser.loadPeaks("chr1");
        assertEquals(14, peaks.size());
    }
}
