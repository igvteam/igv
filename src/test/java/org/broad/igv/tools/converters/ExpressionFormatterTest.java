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
        String outputPath = TestUtils.TMP_OUTPUT_DIR + "testformat.gct";
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
