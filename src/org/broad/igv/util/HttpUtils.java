/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.util;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.util.stream.IGVUrlHelper;
import org.broad.tribble.util.SeekableHTTPStream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URL;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 9/22/11
 */
public abstract class HttpUtils {

    private static Logger log = Logger.getLogger(HttpUtils.class);

    public static boolean byteRangeTested = false;
    public static boolean useByteRange = true;

    public static HttpUtils getInstance() {

        boolean genomeSpaceEnabled = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.GENOME_SPACE_ENABLE);
        return genomeSpaceEnabled ? IGVHttpClientUtils.getInstance() : HttpURLConnectionUtils.getInstance();
    }

    public static boolean isURL(String string) {
        String lcString = string.toLowerCase();
        return lcString.startsWith("http://") || lcString.startsWith("https://") || lcString.startsWith("ftp://")
                || lcString.startsWith("file://");
    }

    /**
     * Test to see if this client can successfully retrieve a portion of a remote file using the byte-range header.
     * This is not a test of the server, but the client.  In some environments the byte-range header gets removed
     * by filters after the request is made by IGV.
     *
     * @return
     */
    public static boolean testByteRange() {

        try {
            String testURL = "http://www.broadinstitute.org/igvdata/byteRangeTest.txt";
            byte[] expectedBytes = {(byte) 'k', (byte) 'l', (byte) 'm', (byte) 'n', (byte) 'o'};

            SeekableHTTPStream str = new SeekableHTTPStream(new IGVUrlHelper(new URL(testURL)));
            str.seek(10);
            byte[] buffer = new byte[5];
            str.read(buffer, 0, 5);

            for (int i = 0; i < buffer.length; i++) {
                if (buffer[i] != expectedBytes[i]) {
                    return false;
                }
            }
            return true;
        } catch (IOException e) {
            log.error("Error while testing byte range ", e);
            // We could not reach the test server, so we can't know if this client can do byte-range tests or
            // not.  Take the "optimistic" view.
            return true;
        }
    }

    public static boolean useByteRange() {
        useByteRange = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.USE_BYTE_RANGE);
        if (useByteRange && !byteRangeTested) {
            useByteRange = testByteRange();
            byteRangeTested = true;
        }
        return useByteRange;
    }

    public abstract void shutdown();

    public abstract String getContentsAsString(URL url) throws IOException;

    public abstract InputStream openConnectionStream(URL url) throws IOException;

    public abstract InputStream openConnectionStream(URL url, boolean abortOnClose) throws IOException;

    public abstract InputStream openConnectionStream(URL url, Map<String, String> headers) throws IOException;

    public abstract InputStream openConnectionStream(URL url, boolean abortOnClose, Map<String, String> headers) throws IOException;

    public abstract boolean resourceAvailable(URL url);

    public abstract String getHeaderField(URL url, String key) throws IOException;

    public abstract long getContentLength(URL url) throws IOException;

    public abstract void updateProxySettings();

   public abstract  boolean downloadFile(String url, File outputFile) throws IOException;


    public abstract void uploadFile(URI uri, File localFile, Map<String,String> headers) throws IOException;;

    public abstract String createDirectory(URL url, String body) throws IOException;;
}
