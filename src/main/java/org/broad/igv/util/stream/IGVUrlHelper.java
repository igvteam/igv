package org.broad.igv.util.stream;

import htsjdk.tribble.util.URLHelper;
import org.apache.log4j.Logger;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

public class IGVUrlHelper implements URLHelper {

    static Logger log = Logger.getLogger(IGVUrlHelperFactory.class);

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
            log.trace("HttpUtils.createURL: " + nurl);
        }
        return HttpUtils.createURL(nurl);
    }

    public boolean exists() {
        log.info("Checking resoure " + url.toExternalForm());
        boolean exists = HttpUtils.getInstance().resourceAvailable(url.toExternalForm());
        log.info("Exists: " + exists);
        return exists;
    }
}
