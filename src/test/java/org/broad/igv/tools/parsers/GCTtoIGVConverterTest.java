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

package org.broad.igv.tools.parsers;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.tools.converters.GCTtoIGVConverter;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GCTtoIGVConverterTest extends AbstractHeadlessTest {

    @Test
    public void testDescriptionMapping() throws IOException {
        String gctFile = TestUtils.DATA_DIR + "gct/igv_test2.gct";
        File igvFile = new File(TestUtils.DATA_DIR + "gct/igv_test2.gct.igv");
        igvFile.deleteOnExit();

        ResourceLocator locator = new ResourceLocator(gctFile);
        GCTtoIGVConverter.convert(locator, igvFile, null, 50000, null, genome);
    }


    @Test
    public void testAffyMapping() throws IOException {
        String gctFile = TestUtils.DATA_DIR + "gct/affy_human.gct";
        File igvFile = new File(TestUtils.DATA_DIR + "gct/affy_human.gct.igv");
        igvFile.deleteOnExit();

        ResourceLocator locator = new ResourceLocator(gctFile);
        GCTtoIGVConverter.convert(locator, igvFile, null, 50000, null, genome);
    }

}
