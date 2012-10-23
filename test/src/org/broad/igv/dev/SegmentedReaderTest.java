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

package org.broad.igv.dev;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * Test of our multi-purpose reader class
 * User: jacob
 * Date: 2012-Oct-23
 */
public class SegmentedReaderTest extends AbstractHeadlessTest {

    @Test
    public void testLoadSegFile() throws Exception {
        String filePath = TestUtils.DATA_DIR + "seg/canFam2_hg18.seg";
        ResourceLocator locator = new ResourceLocator(filePath);
        SegmentedReader reader = new SegmentedReader(locator);
        SegmentedAsciiDataSet ds = reader.loadFromFile();

        assertEquals(18, ds.getChromosomes().size());

    }

    @Ignore
    @Test
    public void testLoadSegDB() throws Exception {
        String filePath = TestUtils.DATA_DIR + "seg/canFam2_hg18.sql";
        ResourceLocator locator = new ResourceLocator(filePath);
        SegmentedReader reader = new SegmentedReader(locator);
        SegmentedAsciiDataSet ds = reader.loadFromDB("canFam2");

        assertEquals(18, ds.getChromosomes().size());

    }


}
