package org.igv.batch;

import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import static org.junit.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 11/30/11
 */
public class TestPortBedgraph {

    private BufferedReader in;

    public TestPortBedgraph() {
    }


    public int importBedGraph(String filename) throws Exception {
        in = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
        String line = "";
        int count = 0;
        while ((line = in.readLine()) != null) {
            count++;
        }
        return count;
    }

    @Test
    public void test1409() throws Exception {
        String testfile = TestUtils.DATA_DIR + "wig/jira_1409.bedgraph";
        int count = importBedGraph(testfile);
        assertTrue(count > 0);
    }


}
