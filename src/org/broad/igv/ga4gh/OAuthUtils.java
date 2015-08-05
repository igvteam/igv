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

package org.broad.igv.ga4gh;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;

/**
 * Created by jrobinso on 11/19/14.
 */
public class OAuthUtils {

    private static Logger log = Logger.getLogger(OAuthUtils.class);

    private static final String REFRESH_TOKEN_KEY = "oauth_refresh_token";
    private static final String PROPERTIES_URL = "https://igvdata.broadinstitute.org/app/oauth_native.json";
    private String genomicsScope = "https://www.googleapis.com/auth/genomics";
    private String gsScope = "https://www.googleapis.com/auth/devstorage.read_only";
    private String profileScope = "https://www.googleapis.com/auth/userinfo.profile";
    private String state = "%2Fprofile";
    private String redirectURI = "http%3A%2F%2Flocalhost%3A60151%2FoauthCallback";
    private String clientId;
    private String clientSecret;
    private String authURI;
    private String tokenURI;

    private String authorizationCode;
    private String accessToken;
    private String refreshToken;
    private long expirationTime;

    public static String GS_HOST = "www.googleapis.com";

    private static OAuthUtils theInstance;
    private String currentUserName;


    public static synchronized OAuthUtils getInstance() {

        if (theInstance == null) {
            theInstance = new OAuthUtils();
        }

        return theInstance;
    }

    private OAuthUtils() {
        // Attempt to fetch refresh token from local store.
        restoreRefreshToken();
    }

    private void fetchOauthProperties() throws IOException {

        String propString = HttpUtils.getInstance().getContentsAsString(new URL(PROPERTIES_URL));
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(propString).getAsJsonObject().get("installed").getAsJsonObject();
        authURI = obj.get("auth_uri").getAsString();
        clientSecret = obj.get("client_secret").getAsString();
        tokenURI = obj.get("token_uri").getAsString();
        clientId = obj.get("client_id").getAsString();
    }

    public void openAuthorizationPage() throws IOException, URISyntaxException {

        if (CommandListener.currentListenerPort != 60151) {
            MessageUtils.showMessage("OAuth failure:  IGV command listener must be open on the default port (60151)");
            return;
        } else {

            if (clientId == null) fetchOauthProperties();

            String url = authURI + "?" +
                    "scope=" + genomicsScope + "%20" + gsScope + "%20" + profileScope + "&" +
                    "state=" + state + "&" +
                    "redirect_uri=" + redirectURI + "&" +
                    "response_type=code&" +
                    "client_id=" + clientId; // Native app

            Desktop.getDesktop().browse(new URI(url));
        }

    }

    // Called from port listener upon receiving the oauth request with a "code" parameter
    public void setAuthorizationCode(String ac) throws IOException {
        authorizationCode = ac;
        fetchTokens();
        fetchUserProfile();
    }

    // Called from port listener upon receiving the oauth request with a "token" parameter TODO -- does this ever happen?
    public void setAccessToken(String accessToken) throws IOException {
        this.accessToken = accessToken;
        fetchUserProfile();
    }

    private void fetchTokens() throws IOException {

        if (clientId == null) fetchOauthProperties();

        URL url = new URL(tokenURI);

        Map<String, String> params = new HashMap<String, String>();
        params.put("code", authorizationCode);
        params.put("client_id", clientId);
        params.put("client_secret", clientSecret);
        params.put("redirect_uri", redirectURI);
        params.put("grant_type", "authorization_code");

        String response = HttpUtils.getInstance().doPost(url, params);
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(response).getAsJsonObject();

        accessToken = obj.getAsJsonPrimitive("access_token").getAsString();
        refreshToken = obj.getAsJsonPrimitive("refresh_token").getAsString();
        expirationTime = System.currentTimeMillis() + (obj.getAsJsonPrimitive("expires_in").getAsInt() * 1000);

        // Try to store in java.util.prefs
        saveRefreshToken();
    }

    /**
     * Fetch a new access token from a refresh token.
     *
     * @throws IOException
     */
    private void refreshAccessToken() throws IOException {

        if (clientId == null) fetchOauthProperties();

        URL url = new URL(tokenURI);

        Map<String, String> params = new HashMap<String, String>();
        params.put("refresh_token", refreshToken);
        params.put("client_id", clientId);
        params.put("client_secret", clientSecret);
        params.put("grant_type", "refresh_token");

        String response = HttpUtils.getInstance().doPost(url, params);
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(response).getAsJsonObject();

        JsonPrimitive atprim = obj.getAsJsonPrimitive("access_token");
        if (atprim != null) {
            accessToken = obj.getAsJsonPrimitive("access_token").getAsString();
            expirationTime = System.currentTimeMillis() + (obj.getAsJsonPrimitive("expires_in").getAsInt() * 1000);
            fetchUserProfile();
        } else {
            // Refresh token has failed, reauthorize from scratch
            reauthorize();
        }

    }

    private void reauthorize() throws IOException {
        logout();
        try {
            openAuthorizationPage();
        } catch (URISyntaxException e) {
            e.printStackTrace();
        }
    }

    private void fetchUserProfile() throws IOException {

        if (clientId == null) fetchOauthProperties();

        URL url = new URL("https://www.googleapis.com/plus/v1/people/me?access_token=" + accessToken);
        String response = HttpUtils.getInstance().getContentsAsJSON(url);
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(response).getAsJsonObject();
        currentUserName = obj.get("displayName").getAsString();
    }


    public String getAccessToken() {

        // Check expiration time, with 1 minute cushion
        if (accessToken == null || (System.currentTimeMillis() > (expirationTime - 60 * 1000))) {
            if (refreshToken != null) {
                try {
                    this.refreshAccessToken();
                } catch (IOException e) {
                    log.error("Error fetching access token", e);
                }
            }
        }

        return accessToken;
    }

    public boolean isLoggedIn() {
        return getAccessToken() != null;
    }

    public String getCurrentUserName() {
        return currentUserName;
    }

    public void logout() {
        accessToken = null;
        refreshToken = null;
        expirationTime = -1;
        currentUserName = null;
        removeRefreshToken();
    }


    private void saveRefreshToken() {
        try {
            Preferences.userRoot().put(REFRESH_TOKEN_KEY, refreshToken);
        } catch (Exception e) {
            log.error("Error storing refresh token", e);
        }
    }


    private void restoreRefreshToken() {
        try {
            refreshToken = Preferences.userRoot().get(REFRESH_TOKEN_KEY, null);
        } catch (Exception e) {
            log.error("Error fetching oauth refresh token", e);
        }
    }


    private void removeRefreshToken() {
        try {
            Preferences.userRoot().remove(REFRESH_TOKEN_KEY);
        } catch (Exception e) {
            log.error("Error removing oauth refresh token", e);
        }
    }


    // Doesn't really belong here....
    public static boolean isGoogleCloud(String url) {
        return url.contains(GS_HOST);
    }

    public void updateSaveOption(boolean aBoolean) {
        if (aBoolean) {
            if (refreshToken != null) {
                saveRefreshToken();
            }
        } else {
            removeRefreshToken();

        }
    }
}
