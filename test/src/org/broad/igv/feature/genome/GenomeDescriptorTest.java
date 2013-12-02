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

package org.broad.igv.feature.genome;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Apr 28, 2011
 */
public class GenomeDescriptorTest {

    /**
     * Test the patch to fix legacy .genome files with incorrect sequence paths
     */

    @Test
    public void testSequencePathFix() {

        // Test for a directory that does exist
        String seqLocation = TestUtils.DATA_DIR;
        GenomeDescriptor desc = new GenomeZipDescriptor("", false, "", "", "", "", "", seqLocation, false, null, null,
                false, false, false, null);
        assertEquals(TestUtils.DATA_DIR, desc.getSequenceLocation());

        // Test for a directory that does not exist
        seqLocation = "/foo/bar";
        desc = new GenomeZipDescriptor("", false, "", "", "", "", "", seqLocation, false, null, null, false, false, false, null);
        assertEquals("/foo/bar", desc.getSequenceLocation());
    }
}
