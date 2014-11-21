package org.broad.igv.ga4gh;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.log4j.Logger;
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
    private static String oauthURL = "https://accounts.google.com/o/oauth2/";
    private static String authEndpoint = "auth";
    private static String tokenEndpoint = "token";
    private static String scope = "https://www.googleapis.com/auth/genomics";
    private static String state = "%2Fprofile";
    private static String clientId = "661332306814-7kotci54n0tr4fdbrff0u79u1m8f7grf.apps.googleusercontent.com";
    private static String clientSecret = "bDph_1LPw3YEZEvHKP2CEBRi";
    private static String redirectURI = "http%3A%2F%2Flocalhost%3A60151%2FoauthCallback";

    private  String authorizationCode;
    private  String accessToken;
    private  String refreshToken;
    private  long expirationTime;

    private static OAuthUtils theInstance;

    public static synchronized OAuthUtils getInstance() {

        if(theInstance == null) {
            theInstance = new OAuthUtils();
        }

        return theInstance;
    }

    private OAuthUtils() {
        // Attempt to fetch refresh token from local store.
        try {
            refreshToken = Preferences.userRoot().get(REFRESH_TOKEN_KEY, null);
        } catch (Exception e) {
            log.error("Error fetching oauth refresh token", e);
        }
    }

    public  void fetchAuthCode() throws IOException, URISyntaxException {

        String url = oauthURL + authEndpoint + "?" +
                "scope=" + scope + "&" +
                "state=" + state + "&" +
                "redirect_uri=" + redirectURI + "&" +
                "response_type=code&" +
                "client_id=" + clientId; // Native app

        Desktop.getDesktop().browse(new URI(url));

    }

    // Called from HttpUtils upon receiving the redirect uri.
    public  void setAuthorizationCode(String ac) throws IOException {
        authorizationCode = ac;
        fetchTokens();
    }

    public  void fetchTokens() throws IOException {

        URL url = new URL(oauthURL + tokenEndpoint);

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
        try {
            Preferences.userRoot().put(REFRESH_TOKEN_KEY, refreshToken);
        } catch (Exception e) {
            log.error("Error storing refresh token", e);
        }
    }

    /**
     * Fetch a new access token from a refresh token.
     *
     * @throws IOException
     */
    public void fetchAccessToken() throws IOException {
        URL url = new URL(oauthURL + tokenEndpoint);

        Map<String, String> params = new HashMap<String, String>();
        params.put("refresh_token", refreshToken);
        params.put("client_id", clientId);
        params.put("client_secret", clientSecret);
        params.put("grant_type", "refresh_token");

        String response = HttpUtils.getInstance().doPost(url, params);
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(response).getAsJsonObject();

        accessToken = obj.getAsJsonPrimitive("access_token").getAsString();
        expirationTime = System.currentTimeMillis() + (obj.getAsJsonPrimitive("expires_in").getAsInt() * 1000);

        // Try to store in java.util.prefs
        try {
            Preferences.userRoot().put(REFRESH_TOKEN_KEY, refreshToken);
        } catch (Exception e) {
            log.error("Error storing refresh token", e);
        }

    }

    public  String getAccessToken() {

        // Check expiration time, with 1 minute cushion
        if(accessToken == null || (System.currentTimeMillis() > (expirationTime - 60*1000))) {
            if(refreshToken != null) {
                try {
                    this.fetchAccessToken();
                } catch (IOException e) {
                    log.error("Error fetching access token", e);
                }
            }
        }

        return accessToken;
    }

    public  void setAccessToken(String accessToken) {
        this.accessToken = accessToken;
    }

    public  boolean isLoggedIn() {
        return getAccessToken() != null;
    }

    public  void logout() {
        accessToken = null;
        refreshToken = null;
        expirationTime = -1;

        try {
            Preferences.userRoot().remove(REFRESH_TOKEN_KEY);
        } catch (Exception e) {
            log.error("Error removing oauth refresh token", e);
        }
    }

    // Main program for testing

    public static void main(String[] args) throws IOException, URISyntaxException {
        getInstance().fetchAuthCode();
    }


}
