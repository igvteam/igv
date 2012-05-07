package org.broad.igv.util;

import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;

import static org.junit.Assert.assertTrue;


/**
 * @author Jim Robinson
 * @date 12/21/11
 */
public class TestByteRange {

    @Test
    public void testByteRange() throws Exception {

        String urlString = "http://www.broadinstitute.org/igv/projects/dev/echo.php";
        String byteRange = "bytes=" + 5 + "-" + 10;
        HttpURLConnection conn = (HttpURLConnection) (new URL(urlString)).openConnection();
        conn.setRequestMethod("GET");
        conn.setRequestProperty("Range", byteRange);

        InputStream is = null;

        int lines = 0;
        try {
            is = conn.getInputStream();
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                lines++;
            }

        } finally {
            if (is != null) {
                is.close();
            }
        }
        assertTrue(lines > 0);
    }
}
