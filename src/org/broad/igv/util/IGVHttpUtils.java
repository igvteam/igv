/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.util;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ftp.FTPClient;
import org.broad.igv.util.ftp.FTPStream;
import org.broad.igv.util.ftp.FTPUtils;
import org.broad.tribble.util.SeekableHTTPStream;

import javax.swing.*;
import java.awt.*;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.*;
import java.util.HashMap;
import java.util.Map;

public class IGVHttpUtils {

    private static Logger log = Logger.getLogger(IGVHttpUtils.class);

    /**
     * Proxy settings (can be null)
     */
    private static ProxySettings proxySettings = null;
    static boolean byteRangeTested = false;
    static boolean useByteRange = true;

    private static boolean testByteRange() {

        try {
            String testURL = "http://www.broadinstitute.org/igvdata/byteRangeTest.txt";
            byte[] expectedBytes = {(byte) 'k', (byte) 'l', (byte) 'm', (byte) 'n', (byte) 'o'};

            SeekableHTTPStream str = new SeekableHTTPStream(new URL(testURL));
            str.seek(10);
            byte[] buffer = new byte[5];
            str.read(buffer, 0, 5);

            for (int i = 0; i < buffer.length; i++) {
                if (buffer[i] != expectedBytes[i]) {
                    return false;
                }
            }
            return true;
        } catch (IOException e) {
            log.error("Error while testing byte range ", e);
            return false;
        }
    }

    public static synchronized boolean useByteRange() {
        useByteRange = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.USE_BYTE_RANGE);
        if (useByteRange && !byteRangeTested) {
            useByteRange = testByteRange();
            byteRangeTested = true;
        }
        return useByteRange;
    }

    // TODO -- show file URLs be included?

    public static boolean isURL(String string) {
        String lcString = string.toLowerCase();
        return lcString.startsWith("http://") || lcString.startsWith("https://") || lcString.startsWith("ftp://")
                || lcString.startsWith("file://");
    }

    public static void disconnect(URLConnection serverConnection) {
        if (serverConnection != null && serverConnection instanceof HttpURLConnection) {
            ((HttpURLConnection) serverConnection).disconnect();
        }
    }


    public static class ProxySettings {
        boolean auth = false;
        String user;
        String pw;
        boolean useProxy;
        String proxyHost;
        int proxyPort = -1;

        public ProxySettings(boolean useProxy, String user, String pw, boolean auth, String proxyHost, int proxyPort) {
            this.auth = auth;
            this.proxyHost = proxyHost;
            this.proxyPort = proxyPort;
            this.pw = pw;
            this.useProxy = useProxy;
            this.user = user;
        }
    }


    /**
     * Wraps url.openConnection(),  and adds proxy authentication if required.
     *
     * @param url
     * @return
     * @throws java.io.IOException
     */


    public static InputStream openConnectionStream(URL url) throws IOException {
        if (url.getProtocol().toLowerCase().equals("ftp")) {
            String userInfo = url.getUserInfo();
            String host = url.getHost();
            String file = url.getPath();
            FTPClient ftp = FTPUtils.connect(host, userInfo);
            ftp.pasv();
            ftp.retr(file);
            return new FTPStream(ftp);
        } else {
            return openHttpStream(url, (Map<String, String>) null);
        }
    }


    public static InputStream openHttpStream(URL url, Map<String, String> requestProperties) throws IOException {

        if(requestProperties == null) {
            requestProperties = new HashMap();
        }
        requestProperties.put("User-Agent", Globals.applicationString());
        HttpURLConnection conn = openHttpConnectionPrivate(url, requestProperties);
        return openHttpStream(url, conn);


    }

    public static InputStream openHttpStream(URL url, HttpURLConnection conn) throws IOException {
        // IF this is a protected directory we will get a 401.  Continue requesting a user / password until
        // the user cancels or the connectino succeeds.

        while (true) {
            InputStream is = null;
            try {
                is = conn.getInputStream();
                return new URLInputStream(conn, is);
            }
            catch (SocketTimeoutException e) {
                throw e;
            }
            catch (IOException e) {
                if (conn.getResponseCode() == 401) {
                    if (is != null) {
                        is.close();
                    }
                    disconnect(conn);
                }
                else {
                    throw e;
                }
            }

            if (getUserPass(url.toExternalForm()) == false) {
                 throw new RuntimeException("A password is required to access " + url.toString());
            }

            Map<String, String> requestProperties = new HashMap();
            for (Map.Entry<String, java.util.List<String>> entry : conn.getRequestProperties().entrySet()) {
                if (entry.getValue().size() > 0) {
                    requestProperties.put(entry.getKey(), entry.getValue().get(0));
                }
            }
            conn.getRequestProperties();
            conn = openHttpConnectionPrivate(url, null);

        }

    }


    private static HttpURLConnection openHttpConnectionPrivate(URL url, Map<String, String> requestProperties) throws IOException {
        URLConnection conn = openConnection(url);
        conn.setConnectTimeout(10000);
        conn.setReadTimeout(60000);

        if (conn instanceof HttpURLConnection) {
            ((HttpURLConnection) conn).setRequestMethod("GET");
        }

        conn.setRequestProperty("Connection", "close");
        if (requestProperties != null) {
            for (Map.Entry<String, String> prop : requestProperties.entrySet()) {
                conn.setRequestProperty(prop.getKey(), prop.getValue());
            }
        }
        return (HttpURLConnection) conn;
    }


    public static URLConnection openConnection(URL url) throws IOException {
        if (useProxy()) {
            Proxy proxy = getProxy();
            URLConnection conn = url.openConnection(proxy);
            if (proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                String encodedUserPwd = base64Encode(proxySettings.user + ":" + proxySettings.pw);
                conn.setRequestProperty("Proxy-Authorization", "Basic " + encodedUserPwd);
                conn.setReadTimeout(60000);
            }
            return conn;
        } else {
            URLConnection conn = url.openConnection();
            return conn;
        }
    }

    public static Proxy getProxy() {
        Proxy proxy = new Proxy(Proxy.Type.HTTP, new InetSocketAddress(proxySettings.proxyHost, proxySettings.proxyPort));
        return proxy;
    }

    public static boolean useProxy() {
        return proxySettings != null && proxySettings.useProxy && proxySettings.proxyHost != null && proxySettings.proxyPort > 0;
    }


    /**
     * Todo,  add comment -- what is returned
     *
     * @param locationString
     * @return
     */
    public static boolean getUserPass(String locationString) {

        //http://www.broadinstitute.org/igvdata/private/cpgIslands.hg18.bed
        JPanel passPanel = new JPanel();
        passPanel.setLayout(new GridLayout(6, 1));

        JLabel message = new JLabel("Please enter your username and password");
        JLabel location = new JLabel(locationString);
        JLabel username = new JLabel("User:");
        JLabel password = new JLabel("Pass:");
        JTextField userField = new JTextField();
        JPasswordField passwordField = new JPasswordField();
        passPanel.add(message);
        passPanel.add(location);
        passPanel.add(username);
        passPanel.add(userField);
        passPanel.add(password);
        passPanel.add(passwordField);

        int a = JOptionPane.showConfirmDialog(IGV.getMainFrame(), passPanel, "Authentication Required", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);

        if (a == JOptionPane.CANCEL_OPTION) {
            return false;
        } else {
            final String userString = userField.getText();
            final char[] userPass = passwordField.getPassword();
            Authenticator.setDefault(new Authenticator() {
                protected PasswordAuthentication getPasswordAuthentication() {
                    return new PasswordAuthentication(userString, userPass);
                }
            });
            return true;
        }
    }

    public static void updateProxySettings() {


        boolean useProxy;
        String proxyHost;
        int proxyPort = -1;
        boolean auth = false;
        String user = null;
        String pw = null;

        PreferenceManager prefMgr = PreferenceManager.getInstance();
        useProxy = prefMgr.getAsBoolean(PreferenceManager.USE_PROXY);
        proxyHost = prefMgr.get(PreferenceManager.PROXY_HOST, null);
        try {
            proxyPort = Integer.parseInt(prefMgr.get(PreferenceManager.PROXY_PORT, "-1"));
        } catch (NumberFormatException e) {
            proxyPort = -1;
        }
        auth = prefMgr.getAsBoolean(PreferenceManager.PROXY_AUTHENTICATE);
        user = prefMgr.get(PreferenceManager.PROXY_USER, null);
        String pwString = prefMgr.get(PreferenceManager.PROXY_PW, null);
        if (pwString != null) {
            pw = Utilities.base64Decode(pwString);
        }


        ProxySettings proxySettings = new ProxySettings(useProxy, user, pw, auth, proxyHost, proxyPort);
        setProxySettings(proxySettings);
    }


    // TODO -- replace use of sun package

    public static String base64Encode(String str) {
        sun.misc.BASE64Encoder encoder = new sun.misc.BASE64Encoder();
        byte[] bytes = str.getBytes();
        return encoder.encode(bytes);

    }

    // TODO -- replace use of sun package

    public static String base64Decode(String str) {
        try {
            return new String((new sun.misc.BASE64Decoder()).decodeBuffer(str));
        } catch (IOException e) {
            log.error("Error decoding string: " + str, e);
            return str;
        }
    }


    public static void setProxySettings(ProxySettings ps) {
        proxySettings = ps;
    }

    public static String getETag(URL url) {
        return getHeaderField(url, "ETag");
    }

    public static String getHeaderField(URL url, String name) {

        // Open an input stream just to check permissions
        InputStream is = null;
        try {
            is = openConnectionStream(url);
        }
        catch (IOException e) {
            log.error("Error getting header field", e);
            return null;
        }
        finally {
            try {
                if (is != null) {
                    is.close();
                }
            } catch (IOException e) {
                log.error("", e);
            }

        }

        URLConnection conn = null;
        try {
            // Create a URLConnection object for a URL
            conn = openConnection(url);
            conn.setReadTimeout(5000);
            return conn.getHeaderField(name);

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (conn != null && conn instanceof HttpURLConnection) {
                ((HttpURLConnection) conn).disconnect();
            }
        }
    }

    public static boolean resourceAvailable(URL url) {
        URLConnection conn = null;

        if (url.getProtocol().toLowerCase().equals("ftp")) {
            return FTPUtils.resourceAvailable(url);
        }

        try {
            // Create a URLConnection object for a URL
            conn = openConnection(url);
            ((HttpURLConnection) conn).setRequestMethod("HEAD");
            conn.setReadTimeout(5000);
            return conn.getHeaderField("ETag") != null;
        } catch (Exception e) {
            return false;
        }
        finally {
            if (conn != null && conn instanceof HttpURLConnection) {
                ((HttpURLConnection) conn).disconnect();
            }
        }
    }

    public static class URLInputStream extends FilterInputStream {
        HttpURLConnection connection;

        protected URLInputStream(HttpURLConnection connection, InputStream inputStream) {
            super(inputStream);
            this.connection = connection;
        }

        @Override
        public void close() throws IOException {
            super.close();
            connection.disconnect();
        }
    }


}
