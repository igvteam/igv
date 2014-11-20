package org.broad.igv.ga4gh;

import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;

/**
 * Created by jrobinso on 11/19/14.
 */
public class GoogleUtils {


    private static String authorizationCode;

    public static void main(String [] args) throws IOException, URISyntaxException {

        String OAUTHURL = "https://accounts.google.com/o/oauth2/auth?";
        String REDIRECT = "http://localhost/igv-web/emptyPage.html";
        String _url = OAUTHURL +
                "scope=https://www.googleapis.com/auth/genomics&" +
                "state=%2Fprofile&" +
                "redirect_uri=http%3A%2F%2Flocalhost%3A60151&" +
                "response_type=code&" +
                "client_id=661332306814-7kotci54n0tr4fdbrff0u79u1m8f7grf.apps.googleusercontent.com"; // Native app
                //"client_id=661332306814-fe8ndvsq0roq3k7e0tvh1bphcms6drvp.apps.googleusercontent.com";  // Web app

       //String response = HttpUtils.getInstance().getContentsAsString(new URL(_url));

        Desktop.getDesktop().browse(new URI(_url));
       // System.out.println(response);

    }

    public static void setAuthorizationCode(String authorizationCode) {
        GoogleUtils.authorizationCode = authorizationCode;
    }

    public static String getAuthorizationCode() {
        return authorizationCode;
    }
}
