package org.broad.igv.util;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;


/**
 * @author Jim Robinson
 * @date 12/21/11
 */
public class TestByteRange {

    public static void main(String[] args) throws IOException {

        String urlString = "http://www.broadinstitute.org/igv/projects/dev/echo.php";
        String byteRange = "bytes=" + 5 + "-" + 10;
        HttpURLConnection conn = (HttpURLConnection) (new URL(urlString)).openConnection();
        conn.setRequestMethod("GET");
        conn.setRequestProperty("Range", byteRange);

        InputStream is = null;
        try {
            is = conn.getInputStream();
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                System.out.println(nextLine);
            }
        } finally {
            if (is != null) {
                is.close();
            }
        }
    }
}
