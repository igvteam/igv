package org.broad.igv.htsget;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Some utilities for interacting with htsget servers
 */
public class HtsgetUtils {

public static Metadata getMetadata(final String url) throws IOException {

            // If the htsget scheme is used, assume the protocol is https, if that fails try http
            if (url.startsWith("htsget://")) {
                try {
                    return getMetadata(url.replace("htsget://", "https://"));
                } catch (IOException e) {
                    return getMetadata(url.replace("htsget://", "http://"));
                }
            }

            URL headerURL = addQueryString(url, "class=header");
            String ticket = HttpUtils.getInstance().getContentsAsJSON(headerURL);
            org.json.JSONObject json;
            try {
                json = new org.json.JSONObject(ticket);
            } catch (org.json.JSONException e) {
                return null;  // Not an htsget server
            }

            if (!json.has("htsget")) {
                // TODO -- throw exception?
                return null;
            } else {
                String format = json.getJSONObject("htsget").getString("format");
                // if(!(format.toUpperCase().equals("BAM")) || format.toUpperCase().equals("VCF")) {
                //     throw new RuntimeException(("Format"))
                // }
                return new Metadata(url, format);
            }
        }

    static URL addQueryString(String urlBase, String queryString) throws MalformedURLException {
        String separator = urlBase.contains("?") ? "&" : "?";
        return new URL(urlBase + separator + queryString);
    }

    public static class Metadata {
        private final String url;
        private final String format;

        public Metadata(String url, String format) {
            this.url = url;
            this.format = format;
        }

        public String getUrl() {
            return url;
        }

        public String getFormat() {
            return format;
        }
    }
}
