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
