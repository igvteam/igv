package org.broad.igv.ga4gh;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;

import static org.junit.Assert.*;

public class Ga4ghTest {



    @Test
    public void testParseReads() throws Exception {

        String testFile = TestUtils.DATA_DIR + "ga4gh/reads.txt";

        BufferedReader br;
        StringBuffer sb = new StringBuffer();
        String line;

        br = new BufferedReader(new FileReader(testFile));
        while((line = br.readLine()) != null) {
            sb.append(line + "\n");
        }

        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(sb.toString()).  getAsJsonObject();

        JsonArray reads = obj.getAsJsonArray("reads");

        System.out.println(obj);

    }

}