/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.tools.parsers;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.converters.GCTtoIGVConverter;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.broad.igv.util.TestUtils.loadGenome;

/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GCTtoIGVConverterTest {

    Genome genome;

    @Before
    public void setUp() throws IOException {
        TestUtils.setUpHeadless();
        genome = loadGenome();
    }

    @After
    public void tearDown() throws Exception {
    }


    @Test
    public void testDescriptionMapping() throws IOException {
        String gctFile = TestUtils.DATA_DIR + "/gct/igv_test2.gct";
        File igvFile = new File(TestUtils.DATA_DIR + "/gct/igv_test2.gct.igv");
        igvFile.deleteOnExit();

        ResourceLocator locator = new ResourceLocator(gctFile);
        GCTtoIGVConverter.convert(locator, igvFile, null, 50000, null, genome);
    }


    @Test
    public void testAffyMapping() throws IOException {
        String gctFile = TestUtils.DATA_DIR + "/gct/affy_human.gct";
        File igvFile = new File(TestUtils.DATA_DIR + "/gct/affy_human.gct.igv");
        igvFile.deleteOnExit();

        ResourceLocator locator = new ResourceLocator(gctFile);
        GCTtoIGVConverter.convert(locator, igvFile, null, 50000, null, genome);
    }

}
