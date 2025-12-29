package org.broad.igv.tdf;

import junit.framework.TestCase;
import org.broad.igv.util.TestUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * @author Jim Robinson
 * @date 5/22/12
 */
public class TDFUtilsTest extends TestCase {


    public void testTdfToBedgraph() throws Exception {

        // Create a test bedgraph file

        File testFile = new File(TestUtils.DATA_DIR + "wig/test.bedgraph");
        String trackLine = "track name=test";

        PrintWriter pw = new PrintWriter(new FileWriter(testFile));
        pw.println(trackLine);
        String chr = "chr1";
        int step = 1000;
        int span = 100;
        int start = 0;
        int end = start + step;

        int nPoints = 10000;
        float[] expectedData = new float[nPoints];
        for (int i = 0; i < nPoints; i++) {
            float expected = (float) Math.random();
            expectedData[i] = expected;
            pw.println(chr + "\t" + start + "\t" + end + "\t" + expected);
        }
        pw.close();

        // Create TDF file


        // Convert to bedgraph


        // Parse and compare results

    }
}
