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

package org.broad.igv.google;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.logging.*;
import org.broad.igv.DirectoryManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.FileUtils;


import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Created by jrobinso on 11/19/14.
 */
public class OAuthUtils {

    private static Logger log = LogManager.getLogger(OAuthUtils.class);

    private static final String PROPERTIES_URL = "https://s3.amazonaws.com/igv.org.app/desktop_google";

    private static OAuthUtils theInstance;

    /**
     * Map of key name -> oauth provider.  This is currently only used to distinguish Amazon from Google providers.
     */
    static Map<String, OAuthProvider> providers;

    static OAuthProvider defaultProvider;

    public static synchronized OAuthUtils getInstance() {
        if (theInstance == null) {
            theInstance = new OAuthUtils();
        }
        return theInstance;
    }


    public OAuthProvider getProvider(String providerName) {
        if (providerName != null) {
            if (!providers.containsKey(providerName)) {
                throw new RuntimeException("Unknwon oAuth provider name: " + providerName);
            }
            return providers.get(providerName);
        }
        return defaultProvider;
    }

    public OAuthProvider getProvider() {
        return defaultProvider;
    }

    private OAuthUtils() {
        providers = new HashMap<>();
        try {
            fetchOauthProperties();
        } catch (Exception e) {
            log.error("Error fetching oAuth properties", e);
        }
    }

    private void fetchOauthProperties() throws IOException {

        // Load a provider config specified in preferences
        String provisioningURL = PreferencesManager.getPreferences().getProvisioningURL();
        log.debug("The provisioning URL from prefs.properties is: " + provisioningURL);
        if (provisioningURL != null && provisioningURL.length() > 0) {
            String json = loadAsString(provisioningURL);
            parseProviderJson(json, provisioningURL);
        }

        // Local config takes precendence, overriding URL provisioned and Broad's default oauth-config.json.gz
        String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config-custom.json";
        if ((new File(oauthConfig)).exists()) {
            try {
                log.debug("Loading Oauth properties from: " + oauthConfig);
                String json = loadAsString(oauthConfig);
                parseProviderJson(json, oauthConfig);
            } catch (IOException e) {
                log.error(e);
            }
        }

        if (defaultProvider == null) {
            // IGV default
            log.debug("$HOME/igv/oauth-config.json not found, reading Java .properties instead from Broad's IGV properties endpoint: " + PROPERTIES_URL);
            String propString = loadAsString(PROPERTIES_URL);
            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(propString).getAsJsonObject();
            defaultProvider = new OAuthProvider(obj);
        }
    }

    public void updateOauthProvider(String provisioningURL) throws IOException {
        if (provisioningURL == null || provisioningURL.trim().length() == 0) {
            defaultProvider = null;
            // TODO -- reset (load google default)?
        } else {
            String json = loadAsString(provisioningURL);
            parseProviderJson(json, provisioningURL);
        }
    }


    private String loadAsString(String urlOrPath) throws IOException {
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
     * @param path
     * @throws IOException
     */
    private void parseProviderJson(String json, String path) throws IOException {
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(json).getAsJsonObject();
        OAuthProvider p = new OAuthProvider(obj);

        // Hack - a move towards multiple providers
        if (obj.has("auth_provider") && obj.get("auth_provider").getAsString().equals("Amazon")) {
            providers.put("Amazon", p);
            AmazonUtils.setCognitoConfig(obj);
        } else {
            if (defaultProvider != null) {
                log.info("Overriding default oAuth provider with " + path);
            }
            defaultProvider = p;
        }
    }

    // Set the authorization code from the callback (redirect URI).  We are sharing a single redirect URI between
    // all providers.  At the moment his is hacked specifically for Google and Amazon, if its not Google we assume
    // its Amazon.  In the future we'll most likely have to use distinct redirect URIs, which means opening more ports.
    public void setAuthorizationCode(Map<String, String> params) throws IOException {
        OAuthProvider provider = null;
        if (params.containsKey("scope") && params.get("scope").contains("googleapis")) {
            provider = defaultProvider;
        } else if (providers.containsKey("Amazon")) {
            provider = providers.get("Amazon");
        }
        if (provider == null) {
            provider = defaultProvider;
        }
        provider.setAuthorizationCode(params.get("code"));
    }


    // Used by batch commands (CommandExecutor).  No argument provided to select provider
    public void setAccessToken(Map<String, String> params) {
        OAuthProvider provider = defaultProvider;
        provider.setAccessToken("token");
    }

    /**
     * Open an input stream for reading a local or remote file.
     * <p>
     * NOTE:  This method does not use HttpUtils by design.  This class must be independent of HttpUtils.
     *
     * @param urlOrPath -- either an http URL or path to a local file.  Can be gzipped
     * @return
     * @throws IOException
     */
    public static InputStream openInputStream(String urlOrPath) throws IOException {

        InputStream inputStream = null;
        if (urlOrPath.startsWith("http://") || urlOrPath.startsWith("https://")) {
            URL url = new URL(urlOrPath);
            URLConnection conn = url.openConnection();
            inputStream = conn.getInputStream();
            urlOrPath = url.getPath();
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
