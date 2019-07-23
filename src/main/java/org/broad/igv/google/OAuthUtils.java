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
import com.google.gson.JsonPrimitive;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.JWTParser;
import software.amazon.awssdk.services.cognitoidentity.model.Credentials;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLDecoder;
import java.time.Duration;
import java.util.*;
import java.util.prefs.Preferences;

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

    private String authProvider = "";
    private String appIdURI = null;
    public static String findString = null;
    public static String replaceString = null;


    private static final String REFRESH_TOKEN_KEY = "oauth_refresh_token";
    private static final String PROPERTIES_URL = "https://s3.amazonaws.com/igv.org.app/desktop_google";
    private String gsScope = "https://www.googleapis.com/auth/devstorage.read_only";
    private String driveScope = "https://www.googleapis.com/auth/drive.readonly";
    private String emailScope = "https://www.googleapis.com/auth/userinfo.email";
    private String state = UUID.randomUUID().toString(); // "RFC6749: An opaque value used by the client to maintain state"
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

    private static OAuthUtils theInstance;
    private String currentUserName;
    private String currentUserEmail;
    private String currentUserID;

    private static JsonObject response;

    // by default this is the google scope
    private String scope = driveScope + "%20" + gsScope + "%20" + emailScope;

    public static synchronized OAuthUtils getInstance() throws IOException {

        if (theInstance == null) {
            theInstance = new OAuthUtils();
        }

        return theInstance;
    }

    private OAuthUtils() throws IOException {
        // XXX: Refactor/rethink this for multiple providers
        //restoreRefreshToken();
        fetchOauthProperties();
    }

    public void fetchOauthProperties() throws IOException {

        String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config.json";

        if (oauthConfig == null || !FileUtils.resourceExists(oauthConfig)) {
            // Remote provisioning of OAUTH attributes
            log.debug("$HOME/igv/oauth-config.json not found, reading Java .properties instead from: "+PROPERTIES_URL);
            String propString = HttpUtils.getInstance().getContentsAsGzippedString(HttpUtils.createURL(PROPERTIES_URL));
            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(propString).getAsJsonObject().get("installed").getAsJsonObject();
            authURI = obj.get("auth_uri").getAsString();
            clientSecret = obj.get("client_secret").getAsString();
            tokenURI = obj.get("token_uri").getAsString();
            clientId = obj.get("client_id").getAsString();
        } else {
            // Local IGV user directory override of OAUTH config
            log.debug("Loading Oauth properties from: " + oauthConfig);
            JsonParser parser = new JsonParser();
            String json = FileUtils.getContents(oauthConfig);
            JsonObject obj = parser.parse(json).getAsJsonObject();

            // Mandatory attributes, fail hard if not present
            try {
                clientId = obj.get("client_id").getAsString();
                authURI = obj.get("authorization_endpoint").getAsString();
                tokenURI = obj.get("token_endpoint").getAsString();
            } catch (Exception e) {
                throw new IOException(oauthConfig+" is missing crucial attributes such as: client_id, " +
                                                  "authorization_endpoint or token_endpoint");
            }
            
            // XXX: What is this for exactly @jrobinso?
            // JsonElement je = obj.get("find_string");
            // if (je != null) {
            //     findString = je.getAsString();
            // }
            // je = obj.get("replace_string");
            // if (je != null) {
            //     replaceString = je.getAsString();

            // Optional or custom attributes, fail on runtime, depending on identity provider configuration
            clientSecret = obj.has("client_secret") ? obj.get("client_secret").getAsString() : null;
            setAuthProvider(obj.has("auth_provider") ? obj.get("auth_provider").getAsString() : authProvider);
            tokenURI = obj.has("token_endpoint") ? obj.get("token_endpoint").getAsString() : null;
            appIdURI = obj.has("app_id_uri") ? obj.get("app_id_uri").getAsString() : null;
            setScope(obj.has("scope") ? obj.get("scope").getAsString() : scope);
            // Accomodate Jim's custom case
            if (scope.equals("none")) {
                setScope(null);
            }

            // @igvteam custom oauth attributes
            GoogleUtils.GS_HOST = obj.has("hosts") ? obj.get("hosts").getAsString() : null;
            authProvider = obj.has("auth_provider") ? obj.get("auth_provider").getAsString() : null;
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

        String redirect = oobURI;
        // if the listener is active, then set the redirect URI.  dwm08
        if (CommandListener.isListening()) {
            redirect = redirectURI;
        }
        String url;
        if (appIdURI == null) {
            log.debug("appIdURI is null, skipping resource setting");
            url = authURI + "?" +
                    "scope=" + scope + "&" +
                    "state=" + state + "&" +
                    "redirect_uri=" + redirect + "&" +
                    "response_type=code&" +
                    "client_id=" + clientId; // Native app
        }

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
    }

    // Called from port listener upon receiving the oauth request with a "token" parameter
    public void setAccessToken(String accessToken) {
        this.accessToken = accessToken;
    }

    public void setScope(String scope) {
        this.scope = scope;
    }

    private void fetchTokens(String redirect) throws IOException {

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

            accessToken = response.get("access_token").getAsString();
            refreshToken = response.get("refresh_token").getAsString();
            idToken = response.get("id_token").getAsString();

            log.debug("Oauth2 refresh token: "+refreshToken);
            log.debug("Oauth2 token_id: "+idToken);
            log.debug("Oauth2 access token: "+accessToken);
            log.debug("Oauth2 state: "+state);

            refreshToken = response.get("refresh_token").getAsString();
            expirationTime = System.currentTimeMillis() + (response.get("expires_in").getAsInt() * 1000);

            JsonObject payload = JWTParser.getPayload(response.get("id_token").getAsString());

            // Populate this class with user profile attributes
            fetchUserProfile(payload);


            if (authProvider.equals("Amazon")) {
                // Get AWS credentials after getting relevant tokens
                Credentials aws_credentials;
                aws_credentials = AmazonUtils.GetCognitoAWSCredentials();

                // Update S3 client with newly acquired token
                AmazonUtils.updateS3Client(aws_credentials);
            }


            // Notify UI that we are authz'd/authn'd
            if (isLoggedIn()) {
                IGVEventBus.getInstance().post(new AuthStateEvent(true, this.authProvider, this.getCurrentUserName()));
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
    public JsonObject fetchUserProfile(JsonObject jwt_payload) throws IOException {
        try {
            currentUserName = jwt_payload.has("name") ? jwt_payload.get("name").getAsString() : null;
            currentUserEmail = jwt_payload.has("email") ? jwt_payload.get("email").getAsString() : null;
            currentUserID = jwt_payload.has("id") ? jwt_payload.get("id").getAsString() : null;

            return jwt_payload;
        } catch (Throwable exception) {
            log.error(exception);
            return null;
        }
    }

    public String getAccessToken() {

        // Check expiration time, with 1 minute cushion
//        if (accessToken == null || (System.currentTimeMillis() > (expirationTime - 60 * 1000))) {
        if (accessToken == null || true) {
            System.out.println("Yay");
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

    /**
    *
    * Get OAuth credential tokens expiration time (in seconds).
    *
    */
    public static Duration getExpirationTime() {
        return Duration.ofMillis(expirationTime - System.currentTimeMillis());
    }

    public static Date getExpirationDate() {
        return new Date(expirationTime);
    }

    public class AuthStateEvent {
        boolean authenticated;
        String authProvider;
        String userName;
        String email;

        // Assuming that if this event is called, we are indeed autz/authn'd
        public AuthStateEvent(boolean authenticated, String authProvider, String userName) {
            this.authenticated = authenticated;
            this.authProvider = authProvider;
            this.userName = userName;
        }

        public boolean isAuthenticated() {
            return authenticated;
        }

        public String getAuthProvider() {
            return authProvider;
        }

        public String getUserName() {
            return userName;
        }

        public String getEmail() {
            return currentUserEmail;
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
                if (GoogleUtils.isGoogleCloud(url)) {

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
