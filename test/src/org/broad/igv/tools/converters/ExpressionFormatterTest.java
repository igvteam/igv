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

package org.broad.igv.tools.converters;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.Dataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertTrue;


/**
 * User: jacob
 * Date: 2012/03/19
 */
public class ExpressionFormatterTest extends AbstractHeadlessTest {

    @Test
    public void testConvertGCT() throws Exception {
        //String inputPath = TestUtils.DATA_DIR + "gct/igv_test2.gct";
        String inputPath = TestUtils.DATA_DIR + "gct/GBM.methylation__sampled.data.txt";
        String outputPath = TestUtils.DATA_DIR + "out/testformat.gct";
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
