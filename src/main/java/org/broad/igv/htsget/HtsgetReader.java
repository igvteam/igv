package org.broad.igv.htsget;

import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * A low-level class to fetch header and data from htsget servers.  No intrepretation or parsing of the header or
 * data is done here, the class is designed as a helpfer for format specific readers.  Currently this only supports the
 * variant (VCF) format, the htsjdk is used for BAM and CRAM which does not use this helper.
 */
class HtsgetReader {

    private String url;
    private String format;

    /**
     * Factory method to return an htsget reader.
     */
    static HtsgetReader getReader(HtsgetUtils.Metadata metadata) {
        String format = metadata.getFormat();
        return new HtsgetReader(metadata.getUrl(), format);
    }

    HtsgetReader(String url, String format) {
        this.url = url;
        this.format = format.toUpperCase();
    }

    byte[] readHeader() throws IOException {
        URL headerURL = HtsgetUtils.addQueryString(this.url, "class=header&format=" + this.format);
        String ticketString = HttpUtils.getInstance().getContentsAsJSON(headerURL);
        org.json.JSONObject ticket = new org.json.JSONObject(ticketString);
        return loadURLs(ticket);
    }

    byte[] readData(String chr, int start, int end) throws IOException {
        final String queryString = "format=" + this.format + "&referenceName=" + chr + "&start=" + start + "&end=" + end;
        URL queryURL = HtsgetUtils.addQueryString(this.url, queryString);
        String ticketString = HttpUtils.getInstance().getContentsAsJSON(queryURL);
        org.json.JSONObject ticket = new org.json.JSONObject(ticketString);
        return loadURLs(ticket);
    }

    private byte[] loadURLs(org.json.JSONObject ticket) throws IOException {

        org.json.JSONObject container = ticket.getJSONObject("htsget");
        org.json.JSONArray urls = container.getJSONArray("urls");
        byte[] bytes = null;
        for (int i = 0; i < urls.length(); i++) {

            org.json.JSONObject next = urls.getJSONObject(i);
            String urlString = next.getString("url");
            if (urlString.startsWith("data:")) {
                throw new RuntimeException("Data URLs are not currently supported");
            }

            URL url = new URL(urlString);
            Map<String, String> headers = new HashMap<>();
            if (next.has("headers")) {
                org.json.JSONObject headerObj = next.getJSONObject("headers");
                for (String key : headerObj.keySet()) {
                    headers.put(key, headerObj.getString(key));
                }
            }

            byte[] b = HttpUtils.getInstance().getContentsAsBytes(url, headers);
            if (bytes == null) {
                bytes = b;
            } else {
                byte[] newBytes = new byte[bytes.length + b.length];
                System.arraycopy(bytes, 0, newBytes, 0, bytes.length);
                System.arraycopy(b, 0, newBytes, bytes.length, b.length);
                bytes = newBytes;
            }
        }
        return bytes;
    }
}

