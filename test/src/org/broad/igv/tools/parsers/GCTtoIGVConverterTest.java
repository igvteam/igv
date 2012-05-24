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
