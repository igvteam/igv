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
