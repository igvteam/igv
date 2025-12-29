package org.broad.igv.util;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.List;
import java.util.Map;

/**
 * Utility class to dump headers from an http response.  Used to debug JWS update problems.
 *
 * @author jrobinso
 * @date May 28, 2011
 */
public class DumpHttpHeaders {

    public static void main(String[] args) throws IOException {

        for(String url : args) {
            dumpFields(url);
        }
    }

    private static void dumpFields(String url) throws IOException {
        URLConnection conn = (HttpUtils.createURL(url)).openConnection();
        Map<String, List<String>> headerFields = conn.getHeaderFields();
        System.out.println(url);
        for (Map.Entry<String, List<String>> entry : headerFields.entrySet()) {
            System.out.print(entry.getKey() + "  ");
            for (String v : entry.getValue()) {
                System.out.println(v + "  ");
            }
        }
        System.out.println();
    }
}
