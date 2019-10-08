package org.broad.igv.util;

import org.apache.log4j.Logger;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

public class URLUtils {

    private static Logger log = Logger.getLogger(URLUtils.class);

    private static Pattern URLmatcher = Pattern.compile(".{1,8}://.*");

    public static String getPath(String url) throws MalformedURLException {
        URL fauxURL = getFauxUrl(url);
        return fauxURL.getPath();
    }

    public static String getQuery(String url) throws MalformedURLException {
        URL fauxURL = getFauxUrl(url);
        return fauxURL.getQuery();
    }


    public static String getHost(String url) throws MalformedURLException {
        URL fauxURL = getFauxUrl(url);
        return fauxURL.getHost();
    }

    public static Map<String, String> parseQueryString(String query) {
        String[] params = query.split("&");
        Map<String, String> map = new HashMap<String, String>();
        for (String param : params) {
            String[] name_val = param.split("=", 2);
            if (name_val.length == 2) {
                map.put(name_val[0], name_val[1]);
            }
        }
        return map;
    }

    public static String addExtension(String url, String extension) {

        try {
            String path = getPath(url);
            return url.replaceFirst(path, path + extension);
        } catch (MalformedURLException e) {
            log.error(e);
            return url + extension;
        }
    }

    public static String replaceExtension(String url, String extension, String newExtension) {

        String path = url;
        try {
            path = getPath(url);
       } catch (MalformedURLException e) {
            log.error(e);
        }
        int idx = path.lastIndexOf(extension);
        String newPath = path.substring(0, idx) + newExtension + path.substring(idx + extension.length());
        return url.replaceFirst(path, newPath);
    }

    public static String addParameter(String urlString, String parameter) {
        return (urlString.indexOf('?') > 0 ? "&" : "?") + parameter;
    }

    /**
     * Checks if the string is a URL (not necessarily remote, can be any protocol)
     *
     * @param f
     * @return
     */
    public static boolean isURL(String f) {
        return URLmatcher.matcher(f).matches();
    }


    /**
     * Convert unsupported protcols (gs, s3, etc) to http for the purpose of parsing out components.
     *
     * @param url
     * @return
     */

    private static URL getFauxUrl(String url) throws MalformedURLException {
        String fauxURL = url;
        int idx = url.indexOf(":");
        if (idx > 0) {
            String protocol = url.substring(0, idx);
            if (!protocol.startsWith("http")) {
                fauxURL = "http" + url.substring(idx);
            }
        }
        return new URL(fauxURL);
    }


}
