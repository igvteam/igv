package org.igv.oauth;

import org.igv.event.IGVEvent;
import org.igv.logging.*;
import org.igv.batch.CommandListener;
import org.igv.event.IGVEventBus;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.util.MessageUtils;
import org.igv.util.AmazonUtils;
import org.igv.util.GoogleUtils;
import org.igv.util.HttpUtils;
import org.igv.util.JWTParser;
import org.json.JSONArray;
import org.json.JSONObject;
import software.amazon.awssdk.services.sts.model.Credentials;

import java.awt.*;
import java.io.IOException;
import java.net.*;
import java.time.Duration;
import java.util.*;

public class OAuthProvider {

    private static Logger log = LogManager.getLogger(OAuthProvider.class);

    private static int TOKEN_EXPIRE_GRACE_TIME = 1000 * 60; // 1 minute

    private String authProvider = "";
    private String appIdURI;
    public static String findString = null;
    public static String replaceString = null;
    private String state; // "RFC6749: An opaque value used by the client to maintain state"
    private String redirectURI;
    private String clientId;
    private String clientSecret;
    private String authEndpoint;
    private String tokenEndpoint;
    private String accessToken;
    private String refreshToken;
    private String codeChallenge;
    private String codeVerifier;
    private String codeChallengeMethod;
    private long expirationTime; // in milliseconds
    private String scope;
    private String[] hosts;
    private String currentUserName;
    private String currentUserID;
    private String currentUserEmail;
    private JSONObject response;


    public OAuthProvider(JSONObject obj) throws IOException {

        state = UUID.randomUUID().toString(); // "RFC6749: An opaque value used by the client to maintain state"

        // For backward compatibility
        if (obj.has("installed") && !obj.has("client_id")) {
            obj = obj.getJSONObject("installed");
        }

        // Mandatory attributes, fail hard if not present
        if (!(obj.has("client_id") &&
                (obj.has("auth_uri") || obj.has("authorization_endpoint")) &&
                (obj.has("token_uri") || obj.has("token_endpoint")))) {
            throw new RuntimeException("oauthConfig is missing crucial attributes such as: client_id, client_secret,  " +
                    "authorization_endpoint/auth_uri, or token_endpoint/token_uri.");
        }

        clientId = obj.getString("client_id");
        authEndpoint = obj.has("auth_uri") ?
                obj.getString("auth_uri") :
                obj.getString("authorization_endpoint");
        tokenEndpoint = obj.has("token_uri") ?
                obj.getString("token_uri") :
                obj.getString("token_endpoint");
        clientSecret = obj.has("client_secret") ? obj.getString("client_secret") : null;

        // Optional attributes
        if (obj.has("auth_provider")) {
            authProvider = obj.getString("auth_provider");
        }
        if (obj.has("scope")) {
            scope = obj.getString("scope");
        } else if (isGoogle()) {
            String gsScope = "https://www.googleapis.com/auth/devstorage.read_only";
            String emailScope = "https://www.googleapis.com/auth/userinfo.email";
            scope = gsScope + "%20" + emailScope;

        }
        appIdURI = obj.has("app_id_uri") ?    // Microsoft Azure.
                obj.getString("app_id_uri") :
                obj.has("resource") ? obj.getString("resource") : null;

        if (obj.has("hosts")) {
            // hosts element may be an array or a single string - put in hosts array either way
            Object hostsObj = obj.get("hosts");
            if (hostsObj instanceof JSONArray) {
                JSONArray hostsArrJson = (JSONArray) hostsObj;
                hosts = new String[hostsArrJson.length()];
                for (int i = 0; i < hostsArrJson.length(); i++)
                    hosts[i] = hostsArrJson.getString(i);
            } else {
                hosts = new String[1];
                hosts[0] = obj.getString("hosts");
            }
        }

        if (obj.has("redirect_uris")) {
            JSONArray urisArray = obj.getJSONArray("redirect_uris");
            redirectURI = urisArray.getString(0);
        } else if (obj.has("redirect_uri")) {
            redirectURI = obj.getString("redirect_uri");
        } else {
            String portNumber = PreferencesManager.getPreferences().getPortNumber();
            redirectURI = "http://localhost:" + portNumber + "/oauthCallback";
        }

        // Generate PKCE challenge and verifier
        codeVerifier = PKCEUtils.generateCodeVerifier();
        try {
            codeChallenge = PKCEUtils.generateCodeChallange(codeVerifier);
            codeChallengeMethod = "S256";
        } catch (Exception e) {
            codeChallenge = codeVerifier;
            codeChallengeMethod = "plain";
            log.error("Error encoding PKCE challenge", e);
        }

        // Deprecated properties -- for backward compatibility
        findString = obj.has("find_string") ? obj.getString("find_string") : null;
        replaceString = obj.has("replace_string") ? obj.getString("replace_string") : null;
    }

    public String getState() {
        return state;
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

        // If the port listener is not on, try starting it
        if (!CommandListener.isListening()) {
            CommandListener.start();
        }

        // if the listener is not active, prompt the user for the access token
        if (!CommandListener.isListening()) {
            String ac = MessageUtils.showInputDialog(
                    "The IGV port listener is off and is required for OAuth authentication through IGV<br/>.  " +
                            "If you have an access token obtained by other means enter it here.");
            if (ac != null) {
                setAccessToken(ac);
            }
        } else {

            String url = authEndpoint +
                    "?state=" + state +
                    "&redirect_uri=" + URLEncoder.encode(redirectURI, "utf-8") +
                    "&client_id=" + clientId +
                    "&response_type=code" +
                    "&code_challenge=" + codeChallenge +
                    "&code_challenge_method=" + codeChallengeMethod;

            if (scope != null) {
                url += "&scope=" + scope;
            }
            if (appIdURI != null) {
                url += "&resource=" + appIdURI;
            }

            // check if the "browse" Desktop action is supported (many Linux DEs cannot directly launch browsers)
            // otherwise, display a dialog box for the user to copy the URL manually.
            Desktop desktop = Desktop.getDesktop();
            if (desktop.isSupported(Desktop.Action.BROWSE)) {
                desktop.browse(new URI(url));
            } else {
                OAuthURLForm.open(IGV.getInstance().getMainFrame(), url);
            }
        }
    }

    // Called from port listener (org.igv.batch.CommandListener) upon receiving the oauth request with a "code" parameter
    public void fetchAccessToken(String authorizationCode) throws IOException {

        URL url = HttpUtils.createURL(tokenEndpoint);

        Map<String, String> params = new HashMap<>();
        params.put("code", authorizationCode);
        params.put("client_id", clientId);
        if (clientSecret != null) {
            params.put("client_secret", clientSecret);
        }
        params.put("state", state);
        params.put("redirect_uri", redirectURI);
        params.put("grant_type", "authorization_code");
        params.put("code_verifier", codeVerifier);

        //  set the resource if necessary for the auth provider
        if (appIdURI != null) {
            params.put("resource", appIdURI);
        }

        try {

            String res = HttpUtils.getInstance().doPost(url, params);
            response = new JSONObject(res);
            accessToken = response.getString("access_token");
            refreshToken = response.getString("refresh_token");
            expirationTime = System.currentTimeMillis() + response.getInt("expires_in") * 1000;
            // Populate this class with user profile attributes
            if (response.has("id_token")) {
                JSONObject payload = JWTParser.getPayload(response.getString("id_token"));
                fetchUserProfile(payload);
            }

            // Notify UI that we are authz'd/authn'd
            if (isLoggedIn()) {
                IGVEventBus.getInstance().post(new AuthStateEvent(true, this.authProvider, this.getCurrentUserName()));
            }

        } catch (Exception e) {
            log.error(e);
        }
    }

    public void setAccessToken(String accessToken) {
        this.accessToken = accessToken;
    }

    /**
     * Fetch a new access token from a refresh token.  Unlike authorization, this is a synchronous operation
     *
     * @throws IOException
     */
    private void refreshAccessToken() throws IOException {

        log.debug("Refresh access token");

        Map<String, String> params = new HashMap<String, String>();
        params.put("refresh_token", refreshToken);
        params.put("client_id", clientId);
        if (clientSecret != null) {
            params.put("client_secret", clientSecret);
        }
        params.put("grant_type", "refresh_token");

        // set the resource if it necessary for the auth provider dwm08
        if (appIdURI != null) {
            params.put("resource", appIdURI);
        }

        // Poke the token refresh endpoint to get new access key
        URL url = HttpUtils.createURL(tokenEndpoint);

        String responseString = HttpUtils.getInstance().doPost(url, params);
        response = new JSONObject(responseString);

        if (response.has("access_token")) {
            accessToken = response.getString("access_token");
            if (response.has("refresh_token")) {
                refreshToken = response.getString("refresh_token");
            }
            expirationTime = System.currentTimeMillis() + response.getInt("expires_in") * 1000;
        } else {
            // Refresh token has failed, reauthorize from scratch
            logout();
            try {
                openAuthorizationPage();
            } catch (URISyntaxException e) {
                e.printStackTrace();
            }
        }
    }


    /**
     * Extract user information from the claim information
     * dwm08
     *
     * @throws IOException
     */
    public JSONObject fetchUserProfile(JSONObject jwt_payload) {
        try {
            currentUserName = jwt_payload.has("name") ? jwt_payload.getString("name") : null;
            if (currentUserName == null && jwt_payload.has("cognito:username")) {
                currentUserName = jwt_payload.getString("cognito:username");
            }
            currentUserEmail = jwt_payload.has("email") ? jwt_payload.getString("email") : null;
            currentUserID = jwt_payload.has("id") ? jwt_payload.getString("id") : null;

            return jwt_payload;
        } catch (Throwable exception) {
            log.error(exception);
            return null;
        }
    }

    public String getAccessToken() {

        // Check expiration time, with 1 minute cushion
        if (accessToken == null || (System.currentTimeMillis() > (expirationTime - TOKEN_EXPIRE_GRACE_TIME))) {
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

    public boolean isLoggedIn() {
        return accessToken != null;
    }

    public String getCurrentUserName() {
        return currentUserName != null ? currentUserName : (currentUserEmail != null ? currentUserEmail : currentUserID);
    }

    public String getCurrentUserEmail() {
        return currentUserEmail;
    }

    public void logout() {
        accessToken = null;
        refreshToken = null;
        expirationTime = -1;
        currentUserName = null;
        IGVEventBus.getInstance().post(new AuthStateEvent(false, this.authProvider, null));
    }

    public JSONObject getAuthorizationResponse() {

        if (response == null) {
            // Go back to auth flow, not auth'd yet
            checkLogin();
        }
        return response;
    }

    /**
     * If not logged in, attempt to login
     */
    public synchronized void checkLogin() {
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
        // wait until authentication successful or 2 minutes -
        // dwm08
        int i = 0;
        while (!isLoggedIn() && i < 1200) {
            ++i;
            try {
                Thread.sleep(100);
            } catch (InterruptedException e1) {
                e1.printStackTrace();
            }
        }
    }

    /**
     * Does this ouath provider apply (should it's access token be used) for the url provided
     *
     * @param url
     * @return
     */
    public boolean appliesToUrl(URL url) {
        // If this provider has a list of hosts, use them to check the url
        if (this.hosts != null && this.hosts.length > 0) {
            for (String host : hosts) {
                if (url.getHost() != null && url.getHost().equals(host)) {
                    return true;
                }
            }
        }
        if (this.isGoogle()) {
            return GoogleUtils.isGoogleURL(url.toExternalForm());
        }
        return false;
    }

    public boolean isGoogle() {
        return this.authEndpoint.contains("accounts.google.com");
    }

    public void setAuthProvider(String authProvider) {
        this.authProvider = authProvider;
    }

    public String getAuthProvider() {
        return authProvider;
    }

    // Assuming that if this event is called, we are indeed autz/authn'd
    public record AuthStateEvent(boolean authenticated, String authProvider, String userName) implements IGVEvent {
    }
}



