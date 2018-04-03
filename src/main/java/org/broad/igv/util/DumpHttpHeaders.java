/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
        URLConnection conn = (new URL(url)).openConnection();
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
