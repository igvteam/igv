package org.igv.tools.converters;

import org.igv.AbstractHeadlessTest;
import org.igv.data.Dataset;
import org.igv.data.expression.ExpressionFileParser;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
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
