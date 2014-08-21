package org.broad.igv.ga4gh;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.sam.Alignment;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Ga4ghTest {



    @Test
    public void testParseReads() throws Exception {

        String testFile = TestUtils.DATA_DIR + "ga4gh/reads.ga4gh.txt";

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

        Iterator<JsonElement> iter = reads.iterator();
        List<Alignment> alignmentList = new ArrayList<Alignment>();
        while(iter.hasNext()) {
            JsonElement next = iter.next();
            Ga4ghAlignment alignment = new Ga4ghAlignment(next.getAsJsonObject());
            alignmentList.add(alignment);
        }
        System.out.println(obj);

    }

}