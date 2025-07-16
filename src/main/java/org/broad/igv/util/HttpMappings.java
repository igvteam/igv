package org.broad.igv.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import org.broad.igv.logging.*;

public class HttpMappings {

    private static Logger log = LogManager.getLogger(HttpMappings.class);
    private static Map<String, String> mappedURLCache = new HashMap<>();

    static {
        mappedURLCache.put("https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/data/url_mappings.tsv",
                "https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/data/url_mappings.tsv");
    }

    static Map<String, String> urlMappings = new HashMap<>();

    private HttpMappings() {
        // Prevent instantiation
    }


    /**
     * Map a URL, for example from a deprecated host or non-http scheme, to a stable newer form. Sort of a pre-request
     * redirect.  This should be used to map to a stable, long-term URL, not to for example time-limited signed URLs.
     *
     * @param urlString
     * @return
     * @throws MalformedURLException
     */
    public static String mapURL(String urlString) throws MalformedURLException {

        // Check cache to avoid unnecessary lookups
        if (mappedURLCache.containsKey(urlString)) {
            return mappedURLCache.get(urlString);
        }

        if (urlMappings.isEmpty()) {
            loadMappings();
        }

        urlString = checkStaticMappings(urlString);
        String key = urlString.startsWith("s3://") || urlString.contains("amazonaws.com") ? getAmazonKey(urlString) : urlString;

        String mappedURL = urlMappings.containsKey(key) ? urlMappings.get(key) : urlString;
        mappedURLCache.put(urlString, mappedURL);       // Record even if not mapped to prevent further lookups
        return mappedURL;
    }

    private static void loadMappings() {

        try {
            URL fileUrl = new URL("https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/data/url_mappings.tsv");
            try (InputStream inputStream = fileUrl.openStream();
                 BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("#") || line.trim().isEmpty()) continue; // Skip comment lines
                    String[] columns = line.split("\t");
                    if (columns.length == 2) {
                        String key = columns[0].trim();
                        // For Amazon S3 URLs use the item's path.  Important because S3 urls can be in global or regional form.
                        if (key.startsWith("s3://") || key.contains("amazonaws.com")) {
                            key = getAmazonKey(key);
                        }
                        urlMappings.put(key, columns[1].trim());
                    }
                }
            }
        } catch (IOException e) {
            log.error("Error loading URL mappings", e);
        }
    }

    /*
     * Convert an Amazon S3 URL to a key to account for global vs regional forms.
     */
    private static String getAmazonKey(String url) throws MalformedURLException {
        return "S3::: " + (new URL(url.replace("s3://", "https://")).getPath());
    }

    private static String checkStaticMappings(String urlString) throws MalformedURLException {
        if (urlString.startsWith("htsget://")) {
            urlString = urlString.replace("htsget://", "https://");
        } else if (urlString.startsWith("gs://")) {
            urlString = GoogleUtils.translateGoogleCloudURL(urlString);
        }

        if (GoogleUtils.isGoogleURL(urlString)) {
            if (urlString.indexOf("alt=media") < 0) {
                urlString = URLUtils.addParameter(urlString, "alt=media");
            }
        }
        String host = URLUtils.getHost(urlString);
        if (host.equals("igv.broadinstitute.org")) {
            urlString = urlString.replace("igv.broadinstitute.org", "s3.amazonaws.com/igv.broadinstitute.org");
        } else if (host.equals("igvdata.broadinstitute.org")) {
            urlString = urlString.replace("igvdata.broadinstitute.org", "s3.amazonaws.com/igv.broadinstitute.org");
        } else if (host.equals("dn7ywbm9isq8j.cloudfront.net")) {
            urlString = urlString.replace("dn7ywbm9isq8j.cloudfront.net", "s3.amazonaws.com/igv.broadinstitute.org");
        } else if (host.equals("www.broadinstitute.org")) {
            urlString = urlString.replace("www.broadinstitute.org/igvdata", "data.broadinstitute.org/igvdata");
        } else if (host.equals("www.dropbox.com")) {
            urlString = urlString.replace("//www.dropbox.com", "//dl.dropboxusercontent.com");
        } else if (host.equals("drive.google.com")) {
            urlString = GoogleUtils.driveDownloadURL(urlString);
        } else if (host.equals("igv.genepattern.org")) {
            urlString = urlString.replace("//igv.genepattern.org", "//igv-genepattern-org.s3.amazonaws.com");
        }

        // data.broadinstitute.org requires https
        urlString = urlString.replace("http://data.broadinstitute.org", "https://data.broadinstitute.org");

        return urlString;
    }

}
