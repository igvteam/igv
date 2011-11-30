package org.broad.igv.util;

import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.io.IOException;
import java.net.*;

/**
 * @author Jim Robinson
 * @date 11/30/11
 */
public class GSTest {


    public static void main(String[] args) throws Exception {

        Authenticator.setDefault(new TestAuthenticator());

        URL url = new URL("https://dm.genomespace.org/datamanager/v1.0/uploadurl/users/jtr/test_session.xml?Content-Length=894&Content-MD5=MI6Hp%2Fjfu4pNt6i8tBNu0A%3D%3D&Content-Type=application/text");
        HttpURLConnection connection = (HttpURLConnection) url.openConnection();

        String token = "XaU32XBOShYZ92Q9zB7txovYaO1UtOxh";
        connection.setRequestProperty("Cookie", "gs-token=" + token);
        connection.setRequestProperty("Accept", "application/json,text/plain");
        connection.setRequestMethod("GET");

        int code = connection.getResponseCode();

        System.out.println("Code = " + code);

    }

    public static class TestAuthenticator extends Authenticator {

        protected PasswordAuthentication getPasswordAuthentication() {

            final String userString = "jtr";
            final char[] userPass = "abc".toCharArray();

            return new PasswordAuthentication(userString, userPass);
        }
    }
}
