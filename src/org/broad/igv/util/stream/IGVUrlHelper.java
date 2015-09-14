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

package org.broad.igv.util.stream;

import org.apache.log4j.Logger;
import org.broad.igv.util.HttpUtils;
import htsjdk.tribble.util.URLHelper;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Jul 6, 2011
 */
public class IGVUrlHelper implements URLHelper {

    static Logger log = Logger.getLogger(IGVUrlHelper.class);

    URL url;

    public IGVUrlHelper(URL url) {
        this.url = url;
    }

    public URL getUrl() {
        return url;
    }

    //Basic caching
    private static Map<URL, Long> contentLengths = new HashMap<URL, Long>();

    public long getContentLength() throws IOException {
        if (contentLengths.containsKey(url)) {
            return contentLengths.get(url);
        } else {
            long length = HttpUtils.getInstance().getContentLength(url);
            contentLengths.put(url, length);
            return length;
        }
    }

    public InputStream openInputStream() throws IOException {
        return HttpUtils.getInstance().openConnectionStream(url);
    }

    public InputStream openInputStreamForRange(long start, long end) throws IOException {

        String byteRange = "bytes=" + start + "-" + end;
        Map<String, String> params = new HashMap();
        params.put("Range", byteRange);
        //Hack for web services which strip range header
        URL url = addStartEndQueryString(start, end);
        return HttpUtils.getInstance().openConnectionStream(url, params);
    }

    /**
     * Add query parameters which should more properly be in Range header field
     * to query string
     *
     * @param start start byte
     * @param end   end byte
     * @throws MalformedURLException
     */
    private URL addStartEndQueryString(long start, long end) throws MalformedURLException {

        String surl = url.toExternalForm();
        String nurl = surl;

        String toadd = String.format("start=%d&end=%d", start, end);
        String[] parts = surl.split("\\?", 2);
        //TODO For now only mess with the string if we already have query parameters
        //nurl = String.format("%s?%s", parts[0], toadd);
        if (parts.length == 2) {
            nurl = String.format("%s?%s", parts[0], toadd);
            nurl += "&" + parts[1];
        }
        if (log.isTraceEnabled()) {
            log.trace("old url: " + surl);
            log.trace("new url: " + nurl);
        }
        return new URL(nurl);
    }

    public boolean exists() {
        return HttpUtils.getInstance().resourceAvailable(url);
    }
}
