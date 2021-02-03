package org.broad.igv.util.blat;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import static org.junit.Assert.*;

public class BlatClientTest {

    @Test
    public void parseUCSCResult() throws Exception {
        String testPath = TestUtils.DATA_DIR + "blat/UCSC_blat_results.html";
        String response = new String(Files.readAllBytes(Paths.get(testPath)));
        List<String> results = BlatClient.parseResult(response);
        assertEquals(5, results.size());
        //System.out.println(results);
    }

    @Test
    public void parseCustomResult() throws Exception {
        String testPath = TestUtils.DATA_DIR + "blat/CUSTOM_blat_results.html";
        String response = new String(Files.readAllBytes(Paths.get(testPath)));
        List<String> results = BlatClient.parseResult(response);
        assertEquals(8, results.size());
        //System.out.println(results);
    }

//    public static void main(String [] args) throws Exception {
//        (new BlatClientTest()).parseUCSCResult();
//        (new BlatClientTest()).parseCustomResult();
//    }
}