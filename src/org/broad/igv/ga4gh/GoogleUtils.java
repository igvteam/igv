package org.broad.igv.ga4gh;

import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by jrobinso on 11/19/14.
 */
public class GoogleUtils {


    private static String authorizationCode;

    private static String oauthURL = "https://accounts.google.com/o/oauth2/";
    private static String authEndpoint = "auth";
    private static String tokenEndpoint = "token";
    private static String scope = "https://www.googleapis.com/auth/genomics";
    private static String state = "%2Fprofile";
    private static String clientId = "661332306814-7kotci54n0tr4fdbrff0u79u1m8f7grf.apps.googleusercontent.com";
    private static String clientSecret = "bDph_1LPw3YEZEvHKP2CEBRi";
    private static String redirectURI = "http%3A%2F%2Flocalhost%3A60151";
    private static String accessToken;
    private static String refreshToken;
    private static long expirationTime;

    public static void main(String[] args) throws IOException, URISyntaxException {

        getAuthCode();
    }

    // Essentially a callback
    public static void setAuthorizationCode(String authorizationCode) {

        GoogleUtils.authorizationCode = authorizationCode;

    }

    public static void getAuthCode() throws IOException, URISyntaxException {

        String url = oauthURL + authEndpoint + "?" +
                "scope=" + scope + "&" +
                "state=" + state + "&" +
                "redirect_uri=" + redirectURI + "&" +
                "response_type=code&" +
                "client_id=" + clientId; // Native app

        Desktop.getDesktop().browse(new URI(url));

    }


    public static void getTokens() throws IOException {

        URL url = new URL(oauthURL + tokenEndpoint);

        Map<String, String> params = new HashMap<String, String>();
        params.put("code", authorizationCode);
        params.put("client_id", clientId);
        params.put("client_secret", clientSecret);
        params.put("redirect_uri", "http%3A%2F%2Flocalhost%3A60151");
        params.put("grant_type", "authorization_code");

        String response = HttpUtils.getInstance().doPost(url, params);
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(response).getAsJsonObject();

        accessToken = obj.getAsJsonPrimitive("access_token").getAsString();
        refreshToken = obj.getAsJsonPrimitive("refresh_token").getAsString();
        expirationTime = System.currentTimeMillis() + (obj.getAsJsonPrimitive("expires_in").getAsInt() * 1000);

        System.out.println("access_token = " + accessToken);

    }
}
