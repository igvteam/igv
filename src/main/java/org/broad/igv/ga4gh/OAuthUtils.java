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
import org.broad.igv.DirectoryManager;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.*;
import java.util.prefs.Preferences;

/**
 * Created by jrobinso on 11/19/14.
 */
public class OAuthUtils {

    private static Logger log = Logger.getLogger(OAuthUtils.class);

    private String authProvider = "";
    private String appIdURI = null;
    public static String findString = null;
    public static String replaceString = null;


    private static final String REFRESH_TOKEN_KEY = "oauth_refresh_token";
    private static final String PROPERTIES_URL = "https://s3.amazonaws.com/igv.org.app/desktop_google";
    private String genomicsScope = "https://www.googleapis.com/auth/genomics";
    private String gsScope = "https://www.googleapis.com/auth/devstorage.read_write";
    private String emailScope = "https://www.googleapis.com/auth/userinfo.email";
    private String state = UUID.randomUUID().toString(); // "An opaque value the client adds to the initial request."
    private String redirectURI = "http%3A%2F%2Flocalhost%3A60151%2FoauthCallback";
    private String oobURI = "urn%3Aietf%3Awg%3Aoauth%3A2.0%3Aoob";
    private String clientId;
    private String clientSecret;
    private String authURI;
    private String tokenURI;

    private String authorizationCode;
    private String accessToken;
    private String refreshToken;
    private static long expirationTime;

    public static String GS_HOST = "www.googleapis.com";

    private static OAuthUtils theInstance;
    private String currentUserName;

    private static JsonObject response;

    // by default this is the google scope
    private String scope = genomicsScope + "%20" + gsScope + "%20" + emailScope;

    public static synchronized OAuthUtils getInstance() throws IOException {

        if (theInstance == null) {
            theInstance = new OAuthUtils();
        }

        return theInstance;
    }

    private OAuthUtils() throws IOException {
        restoreRefreshToken();
        fetchOauthProperties();
    }

    public void fetchOauthProperties() throws IOException {

        // XXX: This needs proper cleanup/testing for existence of attributes... bit of a polyfill situation here
        String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config.json";
        //IGVPreferences.getInstance().get(IGVPreferences.OAUTH_CONFIG);

        if (oauthConfig == null || !FileUtils.resourceExists(oauthConfig)) {
            log.debug("$HOME/igv/oauth-config.json not found, reading Java .properties instead from: "+PROPERTIES_URL);
            String propString = HttpUtils.getInstance().getContentsAsGzippedString(HttpUtils.createURL(PROPERTIES_URL));
            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(propString).getAsJsonObject().get("installed").getAsJsonObject();
            authURI = obj.get("auth_uri").getAsString();
            clientSecret = obj.get("client_secret").getAsString();
            tokenURI = obj.get("token_uri").getAsString();
            clientId = obj.get("client_id").getAsString();
        } else {
            // Experimental -- this will change -- dwm08
            log.debug("Loading Oauth properties from: " + oauthConfig);
            JsonParser parser = new JsonParser();
            String json = FileUtils.getContents(oauthConfig);
            JsonObject obj = parser.parse(json).getAsJsonObject();
//            authURI = obj.get("authorization_endpoint").getAsString();
            authURI = obj.get("auth_uri").getAsString();
//            clientSecret = obj.get("client_secret").getAsString();
//            tokenURI = obj.get("token_endpoint").getAsString();
            tokenURI = obj.get("token_uri").getAsString();
            clientId = obj.get("client_id").getAsString();
//            GS_HOST = obj.get("hosts").getAsString();
//            appIdURI = obj.get("app_id_uri").getAsString();
            appIdURI = obj.get("auth_uri").getAsString();
//            authProvider = obj.get("auth_provider").getAsString();
/*            String scope = obj.get("scope").getAsString();
            if (scope.equals("none")) {
                this.scope = null;
            }
            JsonElement je = obj.get("find_string");
            if (je != null) {
                findString = je.getAsString();
            }
            je = obj.get("replace_string");
            if (je != null) {
                replaceString = je.getAsString();
            }
*/
        }
    }

    /**
     * Send request to authorization provider to start the oauth 2.0
     * authorization process. If the listener is up, wait for a callback
     * Otherwise, provide a dialog where user can provide authentication token.
     * This method has been generalized to use any auth provider (originally google only)
     * dwm08
     *
     * @throws IOException
     * @throws URISyntaxException
     */
    public void openAuthorizationPage() throws IOException, URISyntaxException {
        Desktop desktop = Desktop.getDesktop();

        // properties moved to early init dwm08
        //if (clientId == null) fetchOauthProperties();

        String redirect = oobURI;
        // if the listener is active, then set the redirect URI.  dwm08
        if (CommandListener.isListening()) {
            redirect = redirectURI;
        }
        String url;
        // for auth providers that need scope,
        // set the scope parameter
        //if (scope != null) {
        if (appIdURI == null) {
            log.debug("appIdURI is null, skipping resource setting");
            url = authURI + "?" +
                    "scope=" + scope + "&" +
                    "state=" + state + "&" +
                    "redirect_uri=" + redirect + "&" +
                    "response_type=code&" +
                    "client_id=" + clientId; // Native app
        }
        // for auth providers that need the resource setting
        // the the resource parameter
        //else if (appIdURI != null) {
        else {
            log.debug("appIdURI is not null, setting resource= as part of the authURI");
            url = authURI + "?" +
                    "scope=" + scope + "&" +
                    "state=" + state + "&" +
                    "redirect_uri=" + redirect + "&" +
                    "response_type=code&" +
                    "resource=" + appIdURI + "&" +
                    "client_id=" + clientId; // Native app
        }
//        else {
//        	throw new IOException("Either scope or resource must be provided to authenticate.");
//        }

        log.debug("URL for the auth page is: "+url);

        // check if the "browse" Desktop action is supported (many Linux DEs cannot directly
        // launch browsers!)

        if (desktop.isSupported(Desktop.Action.BROWSE)) {
            desktop.browse(new URI(url));
        } else { // otherwise, display a dialog box for the user to copy the URL manually.
            MessageUtils.showMessage("Copy this authorization URL into your web browser: " + url);
        }

        // if the listener is not active, prompt the user
        // for the access token
        if (!CommandListener.isListening()) {
            String ac = MessageUtils.showInputDialog("Please paste authorization code here:");
            if (ac != null) {
                setAuthorizationCode(ac, oobURI);
            }
        }
    }


    // Called from port listener (org.broad.igv.batch.CommandListener) upon receiving the oauth request with a "code" parameter
    public void setAuthorizationCode(String ac) throws IOException {
        setAuthorizationCode(ac, redirectURI);
    }

    public void setAuthorizationCode(String ac, String redirect) throws IOException {
        authorizationCode = ac;
        log.debug("oauth code parameter: "+ac);
        log.debug("url-encoded redirect_uri: "+redirect);
        fetchTokens(redirect);
        fetchUserProfile();
    }

    // Called from port listener upon receiving the oauth request with a "token" parameter TODO -- does this ever happen?
    public void setAccessToken(String accessToken) throws IOException {
        this.accessToken = accessToken;
        fetchUserProfile();
    }

    public void setScope(String scope) throws IOException {
        this.scope = scope;
    }

    private void fetchTokens(String redirect) throws IOException {

        // properties moved to early init dwm08
        //if (clientId == null) fetchOauthProperties();

        URL url = HttpUtils.createURL(tokenURI);

        Map<String, String> params = new HashMap<String, String>();
        params.put("code", authorizationCode);
        params.put("client_id", clientId);
        // pointless to have clientsecret on a publicly distributed client software?
        if (clientSecret != null) { params.put("client_secret", clientSecret); }
        params.put("redirect_uri", new URLDecoder().decode(redirectURI, "utf-8"));
        params.put("grant_type", "authorization_code");

        //  set the resource if it necessary for the auth provider dwm08
        if (appIdURI != null) {
            params.put("resource", appIdURI);
        }


        try {
            String idToken;

            String res = HttpUtils.getInstance().doPost(url, params);
            JsonParser parser = new JsonParser();

            setResponse(parser.parse(res).getAsJsonObject());
            JsonObject response = getResponse();

            accessToken = response.get("access_token").getAsString();
            refreshToken = response.get("refresh_token").getAsString();
            idToken = response.get("id_token").getAsString();

            log.debug("Oauth2 refresh token: "+refreshToken);
            log.debug("Oauth2 token_id: "+idToken);
            log.debug("Oauth2 access token: "+accessToken);
            log.debug("Oauth2 state: "+state);

            refreshToken = response.get("refresh_token").getAsString();
            expirationTime = System.currentTimeMillis() + (response.get("expires_in").getAsInt() * 1000);

            if (authProvider == "Amazon") {
                // Get AWS credentials after getting relevant tokens
                com.amazonaws.services.cognitoidentity.model.Credentials aws_credentials;
                aws_credentials = AmazonUtils.GetCognitoAWSCredentials();

                // Update S3 client with newly acquired token
                AmazonUtils.updateS3Client(aws_credentials);
            }

            // Notify UI that we are authz'd/authn'd
            if (isLoggedIn()) {
                IGVEventBus.getInstance().post(new AuthStateEvent());
            }

        } catch(Exception e) {
            log.error(e);
            e.printStackTrace();
        }

        // Try to store in java.util.prefs
        saveRefreshToken();
    }



    /**
     * Fetch a new access token from a refresh token.
     *
     * @throws IOException
     */
    private void refreshAccessToken() throws IOException {

        // properties moved to early init dwm08
        //if (clientId == null) fetchOauthProperties();

        URL url = HttpUtils.createURL(tokenURI);

        Map<String, String> params = new HashMap<String, String>();
        params.put("refresh_token", refreshToken);
        params.put("client_id", clientId);
        params.put("client_secret", clientSecret);
        params.put("grant_type", "refresh_token");

        // set the resource if it necessary for the auth provider dwm08
        if (appIdURI != null) {
            params.put("resource", appIdURI);
        }

        String response = HttpUtils.getInstance().doPost(url, params);
        JsonParser parser = new JsonParser();

        setResponse(parser.parse(response).getAsJsonObject());
        JsonObject obj = getResponse();

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

    /**
     * Check if the username is in the claim information. If so, extract it.
     * Otherwise call out to the server to get the current user name.
     * dwm08
     *
     * @throws IOException
     */
    public JsonObject fetchUserProfile() throws IOException {

        try {
            // XXX: This shouldn't belong here, way too Google-specific/hardcoded
            URL url = new URL("https://www.googleapis.com/oauth2/v1/userinfo?access_token=" + accessToken);
            String response = HttpUtils.getInstance().getContentsAsJSON(url);
            JsonParser parser = new JsonParser();

            JsonObject json = parser.parse(response).getAsJsonObject();

            currentUserName = json.get("name").getAsString();
            //currentUserEmail = json.get("email").getAsString();
            //currentUserID = json.get("id").getAsString();

            return json;
        } catch (Throwable exception) {
            log.error(exception);
            return null;
        }
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

    public static long getExpirationTime() {
        return expirationTime;
    }

    public static class AuthStateEvent {
        boolean authenticated;

        // Assuming that if this event is called, we are indeed autz/authn'd
        public AuthStateEvent() {
            this.authenticated = true;
        }

        public boolean isAuthenticated() {
            return this.authenticated;
        }
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
        return url.startsWith("gs://") || url.contains(GS_HOST);
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

    /**
     * Try to login to secure server. dwm08
     */
    public void doSecureLogin() {
        // if user is not currently logged in, attempt to
        // log in user if not logged in dwm08
        if (!isLoggedIn()) {
            try {
                OAuthUtils.getInstance().openAuthorizationPage();
            } catch (Exception ex) {
                MessageUtils.showErrorMessage("Error fetching oAuth tokens.  See log for details", ex);
                log.error("Error fetching oAuth tokens", ex);
            }

        }
        // wait until authentication successful or 1 minute -
        // dwm08
        int i = 0;
        while (!isLoggedIn() && i < 600) {
            ++i;
            try {
                Thread.sleep(100);
            } catch (InterruptedException e1) {
                e1.printStackTrace();
            }
        }
    }

    /**
     * Generate a set of all urls in the session file
     *
     * @param sessionPath
     * @return list of urls
     */
    public static Set<String> findUrlsInSessionFile(String sessionPath) {
        BufferedReader br = null;
        HashSet<String> urlSet = new HashSet<>();
        try {
            br = new BufferedReader(new FileReader(new File(sessionPath)));
            String line;
            while ((line = br.readLine()) != null) {
                int start = line.indexOf("http");
                if (start != -1) {
                    int mid = line.indexOf("://", start);
                    int end = line.indexOf("/", mid + 3);
                    String url = line.substring(start, end);
                    urlSet.add(url);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                }
            }
        }
        return urlSet;
    }

    /**
     * Check if any reference in the session file refers to a server protected
     * by the oauth protocol. If so, check to see if the user is logged in. If
     * user is not logged in, put up login prompt.
     *
     * @param sessionPath
     */
    public void checkServerLogin(String sessionPath) {
        Set<String> urlSet = findUrlsInSessionFile(sessionPath);
        if (urlSet.size() > 0) {
            for (String url : urlSet) {
                if (isGoogleCloud(url)) {

                    doSecureLogin();

                    // user is logged in. Can proceed with the load
                    return;
                }
            }
        }
    }

    public static JsonObject getResponse() {
        return response;
    }

    public static void setResponse(JsonObject res) {
        response = res;
    }

    public String getAuthProvider() {
        return authProvider;
    }

    public void setAuthProvider(String authProvider) {
        this.authProvider = authProvider;
    }
}
