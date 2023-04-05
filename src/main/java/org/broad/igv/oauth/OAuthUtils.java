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

package org.broad.igv.oauth;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.logging.*;
import org.broad.igv.DirectoryManager;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.GoogleUtils;
import org.broad.igv.util.HttpUtils;


import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * Created by jrobinso on 11/19/14.
 */
public class OAuthUtils {

    private static Logger log = LogManager.getLogger(OAuthUtils.class);

    private static final String PROPERTIES_URL = "https://s3.amazonaws.com/igv.org.app/desktop_google";

    private static OAuthUtils theInstance;

    static OAuthProvider awsProvider;

    static OAuthProvider defaultProvider;

    static Map<String, OAuthProvider> providerCache;

    public static synchronized OAuthUtils getInstance() {
        if (theInstance == null) {
            theInstance = new OAuthUtils();
        }
        return theInstance;
    }

    private OAuthUtils() {
        try {
            providerCache = new LinkedHashMap<>();   // Ordered (linked) map is important
            fetchOauthProperties();
        } catch (Exception e) {
            log.error("Error fetching oAuth properties", e);
        }
    }

    /**
     * Called by AWS code only
     */
    public OAuthProvider getAWSProvider() {
        if (awsProvider == null) {
            throw new RuntimeException("AWS Oauth is not configured");
        }
        return awsProvider;
    }

    public OAuthProvider getDefaultProvider() {
        return defaultProvider;
    }

    private void fetchOauthProperties() throws IOException {

        // Load a provider config specified in preferences
        String provisioningURL = PreferencesManager.getPreferences().getProvisioningURL();
        log.debug("The provisioning URL from prefs.properties is: " + provisioningURL);
        if (provisioningURL != null && provisioningURL.length() > 0) {
            String json = loadAsString(provisioningURL);
            parseProviderJson(json);
        }

        // Local config takes precendence, overriding URL provisioned and Broad's default
        String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config.json";
        if (!(new File(oauthConfig)).exists()) {
            oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config-custom.json";
        }
        if ((new File(oauthConfig)).exists()) {
            try {
                log.debug("Loading Oauth properties from: " + oauthConfig);
                String json = loadAsString(oauthConfig);
                parseProviderJson(json);
            } catch (IOException e) {
                log.error(e);
            }
        }

        // Default (Google) provider
        if (defaultProvider == null) {
            String propString = loadAsString(PROPERTIES_URL);
            defaultProvider = parseProviderJson(propString);
            if (defaultProvider.getAuthProvider() == null || defaultProvider.getAuthProvider().isEmpty()) {
                defaultProvider.setAuthProvider("Google");
            }
        }
    }

    /**
     * Called from preferences editor
     *
     * @param provisioningURL
     * @throws IOException
     */
    public void updateOauthProvider(String provisioningURL) throws IOException {
        if(provisioningURL != null && provisioningURL.trim().length() > 0) {
            String json = loadAsString(provisioningURL);
            parseProviderJson(json);
        }
    }


    private String loadAsString(String urlOrPath) throws IOException {
        if(HttpUtils.isRemoteURL(urlOrPath)) {
            urlOrPath = HttpUtils.mapURL(urlOrPath);
        }
        InputStream is = null;
        try {
            is = openInputStream(urlOrPath);
            byte[] bytes = is.readAllBytes();
            if ((bytes[0] == (byte) (GZIPInputStream.GZIP_MAGIC)) && (bytes[1] == (byte) (GZIPInputStream.GZIP_MAGIC >> 8))) {
                bytes = new GZIPInputStream(new ByteArrayInputStream(bytes)).readAllBytes();
            }

            return new String(bytes, "UTF-8");
        } finally {
            if (is != null) is.close();
        }
    }

    /**
     * Parse Oauth provider configuration and update state
     * <p>
     * TODO -- refactor to remove side effects.
     *
     * @param json
     * @throws IOException
     */
    private OAuthProvider parseProviderJson(String json) throws IOException {
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(json).getAsJsonObject();
        OAuthProvider p = new OAuthProvider(obj);
        providerCache.put(p.getState(), p);
        if (obj.has("auth_provider") && obj.get("auth_provider").getAsString().equals("Amazon")) {
            awsProvider = p;
            AmazonUtils.setCognitoConfig(obj);
        }
        if (p.isGoogle()) {
            defaultProvider = p;
        }
        return p;
    }

    public OAuthProvider getProviderForState(String state) throws IOException {
        if (providerCache.containsKey(state)) {
            return providerCache.get(state);
        } else {
            // This should never happen, perhaps an error should be thrown.
            log.warn("No oAuth provider found for callback");
            return defaultProvider;
        }
    }

    public OAuthProvider getProviderForURL(URL url) {
        for (OAuthProvider provider : providerCache.values()) {
            if (provider.appliesToUrl(url)) {
                return provider;
            }
        }

        // With a properly configured IGV we should never get here
        if (GoogleUtils.isGoogleURL(url.toExternalForm()) && defaultProvider.isGoogle()) {
            return defaultProvider;
        } else {
            return null;
        }
    }

    /**
     * Open an input stream for reading a local or remote file.

     * @param urlOrPath -- either an http URL or path to a local file.  Can be gzipped
     * @return
     * @throws IOException
     */
    public static InputStream openInputStream(String urlOrPath) throws IOException {

        InputStream inputStream = null;
        if (urlOrPath.startsWith("http://") || urlOrPath.startsWith("https://")) {
            URL url = new URL(urlOrPath);
            URLConnection conn = HttpUtils.getInstance().openProxiedConnection(url);
            inputStream = conn.getInputStream();
        } else {
            if (urlOrPath.startsWith("file://")) {
                urlOrPath = urlOrPath.substring(7);
            }
            File file = new File(urlOrPath);
            inputStream = new FileInputStream(file);
        }

        return inputStream;

    }


}
