/*
 * Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NON-INFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.tools.converters;

import org.broad.igv.data.Dataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertTrue;


/**
 * User: jacob
 * Date: 2012/03/19
 */
public class ExpressionFormatterTest {

    Genome genome;

    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();
        genome = TestUtils.loadGenome();

    }

    @After
    public void tearDown() throws Exception {
        genome = null;
    }

    @Test
    public void testConvertGCT() throws Exception {
        //String inputPath = TestUtils.DATA_DIR + "/gct/igv_test2.gct";
        String inputPath = TestUtils.DATA_DIR + "/gct/GBM.methylation__sampled.data.txt";
        String outputPath = TestUtils.DATA_DIR + "/out/testformat.gct";
        File inputFile = new File(inputPath);
        File outputFile = new File(outputPath);

        outputFile.delete();
        outputFile.deleteOnExit();

        ExpressionFormatter formatter = new ExpressionFormatter();
        formatter.convert(inputFile, outputFile, ExpressionFormatter.FileType.GCT);

        ExpressionFileParser parser = new ExpressionFileParser(new ResourceLocator(outputPath), null, genome);
        Dataset ds = parser.createDataset();
        assertTrue(ds.getDataMax() > 0);
        assertTrue(ds.getDataMin() < 0);
        //Should this be true?
        //assertTrue(ds.isLogNormalized());
    }
}
