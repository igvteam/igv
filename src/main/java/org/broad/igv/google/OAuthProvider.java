package org.broad.igv.google;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.JWTParser;
import software.amazon.awssdk.services.sts.model.Credentials;

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

public class OAuthProvider {

    private static Logger log = Logger.getLogger(OAuthProvider.class);

    private String authProvider = "";
    private String appIdURI = null;
    public static String findString = null;
    public static String replaceString = null;


    private static final String REFRESH_TOKEN_KEY = "oauth_refresh_token";

    private String state = UUID.randomUUID().toString(); // "RFC6749: An opaque value used by the client to maintain state"
    private String portNumber = PreferencesManager.getPreferences().getPortNumber();
    private String redirectURI = "http%3A%2F%2Flocalhost%3A"+portNumber+"%2FoauthCallback";
    private String oobURI = "urn%3Aietf%3Awg%3Aoauth%3A2.0%3Aoob";
    private String clientId;
    private String clientSecret;
    private String authURI;
    private String tokenURI;
    private String authorizationCode;
    private String accessToken;
    private String refreshToken;
    private long expirationTime; // in milliseconds
    private String scope;
    private String currentUserName;
    private String currentUserID;
    private String currentUserEmail;
    private JsonObject response;
    private JsonObject config;


    public OAuthProvider(JsonObject obj) throws IOException {

        config = obj;

        // Mandatory attributes, fail hard if not present
        try {
            clientId = obj.get("client_id").getAsString();
            authURI = obj.has("auth_uri") ?
                    obj.get("auth_uri").getAsString() :
                    obj.get("authorization_endpoint").getAsString();
            tokenURI = obj.has("token_uri") ?
                    obj.get("token_uri").getAsString() :
                    obj.get("token_endpoint").getAsString();

        } catch (Exception e) {
            log.error(e);
            throw new IOException("oauthConfig is missing crucial attributes such as: client_id, " +
                    "authorization_endpoint/auth_uri or token_endpoint/token_uri");
        }

        // Optional or custom attributes, fail on runtime, depending on identity provider configuration
        clientSecret = obj.has("client_secret") ? obj.get("client_secret").getAsString() : null;
        setAuthProvider(obj.has("auth_provider") ? obj.get("auth_provider").getAsString() : authProvider);
        appIdURI = obj.has("app_id_uri") ? obj.get("app_id_uri").getAsString() : null;
        if (obj.has("scope")) {
            scope = obj.get("scope").getAsString();
        }

        // Special Google properties
        if (authURI.contains("google")) {
            if (scope == null) {
                String gsScope = "https://www.googleapis.com/auth/devstorage.read_only";
                String driveScope = "https://www.googleapis.com/auth/drive.readonly";
                String emailScope = "https://www.googleapis.com/auth/userinfo.email";
                scope = gsScope + "%20" + emailScope;

            }
            if (authProvider == null) {
                authProvider = "Google";
            }
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
        } else {
            log.debug("appIdURI is not null, setting resource= as part of the authURI");
            url = authURI + "?" +
                    "scope=" + scope + "&" +
                    "state=" + state + "&" +
                    "redirect_uri=" + redirect + "&" +
                    "response_type=code&" +
                    "resource=" + appIdURI + "&" +
                    "client_id=" + clientId; // Native app
        }

        log.debug("URL for the auth page is: " + url);

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
        log.debug("oauth code parameter: " + ac);
        log.debug("url-encoded redirect_uri: " + redirect);
        fetchTokens(redirect);
    }

    // Called from port listener upon receiving the oauth request with a "token" parameter
    public void setAccessToken(String accessToken) {
        this.accessToken = accessToken;
    }

    private void fetchTokens(String redirect) throws IOException {

        URL url = HttpUtils.createURL(tokenURI);

        Map<String, String> params = new HashMap<String, String>();
        params.put("code", authorizationCode);
        params.put("client_id", clientId);
        // pointless to have clientsecret on a publicly distributed client software?
        if (clientSecret != null) {
            params.put("client_secret", clientSecret);
        }
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

            log.debug("Oauth2 refresh token: " + refreshToken);
            log.debug("Oauth2 token_id: " + idToken);
            log.debug("Oauth2 access token: " + accessToken);
            log.debug("Oauth2 state: " + state);

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

        } catch (Exception e) {
            log.error(e);
            e.printStackTrace();
        }
    }


    /**
     * Fetch a new access token from a refresh token.
     *
     * @throws IOException
     */
    private void refreshAccessToken() throws IOException {

        // properties moved to early init dwm08
        //if (clientId == null) fetchOauthProperties();

        Map<String, String> params = new HashMap<String, String>();
        params.put("refresh_token", refreshToken);
        params.put("client_id", clientId);
        params.put("client_secret", clientSecret);
        params.put("grant_type", "refresh_token");

        // set the resource if it necessary for the auth provider dwm08
        if (appIdURI != null) {
            params.put("resource", appIdURI);
        }

        // Poke the token refresh endpoint to get new access key
        URL url = HttpUtils.createURL(tokenURI);

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
            if (currentUserName == null && jwt_payload.has("cognito:username")) {
                currentUserName = jwt_payload.get("cognito:username").getAsString();
            }
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
        if (accessToken == null || (System.currentTimeMillis() > (expirationTime - Globals.TOKEN_EXPIRE_GRACE_TIME))) {
            log.debug("Refreshing access token!");
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
     * Get OAuth credential tokens expiration time (in seconds).
     */
    public Duration getExpirationTime() {
        Duration expiration = Duration.ofMillis(expirationTime - System.currentTimeMillis());
        log.debug("Current expiration time of credentials (and presigned urls is): " + expiration.toSeconds() + " seconds and expirationTime in class is: " + expirationTime);
        return expiration;
    }

    public JsonObject getConfig() {
        return config;
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
                // SECURITY: this is definitely a bad practice, we should NOT persist (refresh) tokens locally, ever.
                //saveRefreshToken();
            }
        } else {
            removeRefreshToken();

        }
    }

    /**
     * Try to login to secure server. dwm08
     */
    public synchronized void doSecureLogin() {
        // if user is not currently logged in, attempt to
        // log in user if not logged in dwm08
        if (!isLoggedIn()) {
            try {
                openAuthorizationPage();
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

    public JsonObject getResponse() {
        return response;
    }

    public void setResponse(JsonObject res) {
        response = res;
    }

    public String getAuthProvider() {
        return authProvider;
    }

    public void setAuthProvider(String authProvider) {
        this.authProvider = authProvider;
    }

}
