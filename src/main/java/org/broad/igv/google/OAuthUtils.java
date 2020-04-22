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
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

// XXX: Both Oauth and JWT classes need a serious refactor/audit. Multi-provider support, concurrency, security, etc...:
// WARNING: This class requires a thorough security audit of Oauth/JWT implementations. I would recommend going through:
// https://www.owasp.org/index.php/JSON_Web_Token_(JWT)_Cheat_Sheet_for_Java
//
// And potentially refactor/substitute this logic with:
//
// https://github.com/auth0/java-jwt.
// https://developers.google.com/api-client-library/java/google-oauth-java-client/

/**
 * Created by jrobinso on 11/19/14.
 */
public class OAuthUtils {

    private static Logger log = Logger.getLogger(OAuthUtils.class);

    public static String findString = null;
    public static String replaceString = null;
    private static final String PROPERTIES_URL = "https://s3.amazonaws.com/igv.org.app/desktop_google";

    private static OAuthUtils theInstance;
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

    public void fetchOauthProperties() throws IOException {
        // Load a provider config specified in preferences
        String provisioningURL = PreferencesManager.getPreferences().getProvisioningURL();
        log.debug("The provisioning URL from prefs.properties is: "+provisioningURL);
        if (provisioningURL != null && provisioningURL.length() > 0) {
            loadProvisioningURL(provisioningURL);
        }

        // Local config takes precendence, overriding URL provisioned and Broad's default oauth-config.json.gz
        String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config.json";
        if ((new File(oauthConfig)).exists()) {
            try {
                log.debug("Loading Oauth properties from: " + oauthConfig);
                String json = FileUtils.getContents(oauthConfig);
                parseProviderJson(json, oauthConfig);
            } catch (IOException e) {
                log.error(e);
            }
        }

        if (defaultProvider == null) {
            // IGV default
            log.debug("$HOME/igv/oauth-config.json not found, reading Java .properties instead from Broad's IGV properties endpoint: " + PROPERTIES_URL);
            String propString = HttpUtils.getInstance().getContentsAsGzippedString(HttpUtils.createURL(PROPERTIES_URL));
            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(propString).getAsJsonObject().get("installed").getAsJsonObject();
            defaultProvider = new OAuthProvider(obj);
        }
    }

    public void loadProvisioningURL(String provisioningURL) throws IOException {
        if (provisioningURL != null && provisioningURL.length() > 0) {
            InputStream is = ParsingUtils.openInputStream(provisioningURL);
            String json = ParsingUtils.readContentsFromStream(is);
            parseProviderJson(json, provisioningURL);
        }
    }

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


    // jtr -- I don't think this is used.  If it is I don't know how to distinguish provider.
    public void setAccessToken(Map<String, String> params) {
        OAuthProvider provider = defaultProvider;
        provider.setAccessToken("token");
    }
}
