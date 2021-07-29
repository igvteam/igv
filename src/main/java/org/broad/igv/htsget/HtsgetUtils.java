package org.broad.igv.htsget;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonSyntaxException;
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
        JsonParser parser = new JsonParser();
        JsonObject json = null;
        try {
            json = parser.parse(ticket).getAsJsonObject();
        } catch (JsonSyntaxException e) {
            return null;  // Not an htsget server
        }

        if (!json.has("htsget")) {
            // TODO -- throw exception?
            return null;
        } else {
            String format = json.get("htsget").getAsJsonObject().get("format").getAsString();
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
