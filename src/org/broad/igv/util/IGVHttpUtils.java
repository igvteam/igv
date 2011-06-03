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


import org.apache.commons.httpclient.Credentials;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.UsernamePasswordCredentials;
import org.apache.commons.httpclient.auth.AuthScope;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ftp.FTPClient;
import org.broad.igv.util.ftp.FTPStream;
import org.broad.igv.util.ftp.FTPUtils;
import org.broad.tribble.util.SeekableHTTPStream;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.net.*;
import java.util.*;

public class IGVHttpUtils {

    private static Logger log = Logger.getLogger(IGVHttpUtils.class);

    /**
     * Proxy settings (can be null)
     */
    private static ProxySettings proxySettings = null;
    static boolean byteRangeTested = false;
    static boolean useByteRange = true;
    private static Map<String, java.util.List<String>> gsCookies = new HashMap();
    public static  Credentials GENOME_SPACE_CREDS = null;


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

    /**
     * Create a file from an input stream.
     *
     * @param url
     * @param outputFile
     * @throws java.io.IOException
     */
    public static boolean downloadFile(URL url, File outputFile) throws IOException {

        int maxTries = 1000;
        int nTries = 0;

        log.info("Downloading " + url + " to " + outputFile.getAbsolutePath());
        int downloaded = 0;
        byte[] buf = new byte[64 * 1024]; // 64K buffer

        OutputStream out = null;
        InputStream is = null;
        try {

            out = new FileOutputStream(outputFile);

            URLConnection connection = openConnection(url);
            int contentLength = connection.getContentLength();

            if (contentLength <= 0) {
                // We don't know the content-length, Try downloading until a loop returns zero bytes
                log.debug("Content-length = " + contentLength);

                is = connection.getInputStream();
                int bytesRead;
                while ((bytesRead = is.read(buf)) != -1) {
                    out.write(buf, 0, bytesRead);
                    downloaded += bytesRead;
                }


            } else { // Content length is known,  keep trying until we get it all.
                while (downloaded < contentLength && nTries < maxTries) {
                    is = connection.getInputStream();
                    int bytesRead;
                    while ((bytesRead = is.read(buf)) != -1) {
                        out.write(buf, 0, bytesRead);
                        downloaded += bytesRead;
                    }

                    if (contentLength > downloaded) {
                        is.close();
                        connection = openConnection(url);
                        connection.setRequestProperty("Range", "bytes=" + downloaded + "-");
                        nTries++;
                        log.info("Restarting download from position: " + downloaded);
                    }


                }
            }

            if (downloaded < contentLength) {
                out.close();
                out = null;
                outputFile.delete();
                MessageUtils.showMessage("Error downloading file: " + outputFile.getAbsoluteFile());
                return false;
            } else {
                log.info("Download complete.  Transferred " + downloaded + " bytes");
                return true;
            }


        }


        catch (Exception e) {
            out.close();
            out = null;
            outputFile.delete();
            MessageUtils.showMessage("<html>Error downloading file: " + outputFile.getAbsoluteFile() +
                    "<br/>" + e.toString());
            return false;

        }
        finally {
            if (is != null) {
                is.close();
            }
            if (out != null) {
                out.flush();
                out.close();
            }
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
        } else if (url.toString().contains("genomespace.org")) {
            return openGSConnectionStream(url);
        } else {
            return openHttpStream(url);
        }
    }

    public static InputStream openGSConnectionStream(URL url) throws IOException {

        HttpClient client = new HttpClient();

        // We know all GS urls will need credentials, so get them now
        if (GENOME_SPACE_CREDS == null) {
            String userpass = getUserPass(url.toString());
            if (userpass == null) {
                throw new RuntimeException("Access denied to: " + url.toString());
            }
            String[] tokens = userpass.split(":");
            GENOME_SPACE_CREDS = new UsernamePasswordCredentials(tokens[0], tokens[1]);
        }


        client.getState().setCredentials(new AuthScope("identitytest.genomespace.org", 8443, AuthScope.ANY_REALM), GENOME_SPACE_CREDS);

        GetMethod getMethod = new GetMethod(url.toExternalForm());
        getMethod.setDoAuthentication(true);
        client.getParams().setParameter("http.protocol.allow-circular-redirects", true);
        int status = client.executeMethod(getMethod);
        if(status == 401) {
            // Try again
            getMethod.releaseConnection();
            GENOME_SPACE_CREDS = null;

            return openGSConnectionStream( url);
        }
        if(status != 200) {
            getMethod.releaseConnection();
            throw new RuntimeException("Error connecting.  Status code = " + status);
        }
        return new HttpClientInputStream(getMethod, getMethod.getResponseBodyAsStream());

    }


    public static InputStream openHttpStream(URL url) throws IOException {

        HttpURLConnection conn = openConnection(url);
        InputStream is = conn.getInputStream();
        return new URLInputStream(conn, is);
    }

    public static HttpURLConnection openConnection(URL url) throws IOException {

        HttpURLConnection conn = null;

        if (useProxy()) {
            Proxy proxy = getProxy();
            conn = (HttpURLConnection) url.openConnection(proxy);
            if (proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                String encodedUserPwd = base64Encode(proxySettings.user + ":" + proxySettings.pw);
                conn.setRequestProperty("Proxy-Authorization", "Basic " + encodedUserPwd);
            }
        } else {
            conn = (HttpURLConnection) url.openConnection();
        }

        if (conn.getResponseCode() == 401) {

            final String userPass = getUserPass(url.toExternalForm());
            if (userPass == null) {
                throw new RuntimeException("A password is required to access " + url.toString());
            }
            Authenticator.setDefault(new Authenticator() {
                protected PasswordAuthentication getPasswordAuthentication() {
                    String[] tokens = userPass.split(":");
                    return new PasswordAuthentication(tokens[0], tokens[1].toCharArray());
                }
            });


            disconnect(conn);
            return openConnection(url);

        } else {
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


    public static String getUserPass(String locationString) {

        //http://www.broadinstitute.org/igvdata/private/cpgIslands.hg18.bed
        JPanel passPanel = new JPanel();
        passPanel.setLayout(new GridLayout(6, 1));

        JLabel message = new JLabel("Please enter your username and password");
        if (locationString.length() > 80) {
            locationString = "..." + locationString.substring((locationString.length() - 80));
        }
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

        int a = JOptionPane.showConfirmDialog(null, passPanel, "Authentication Required", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);

        if (a == JOptionPane.CANCEL_OPTION) {
            return null;
        } else {
            final String userString = userField.getText();
            final String userPass = new String(passwordField.getPassword());
            return userString + ":" + userPass;
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

    public static long getContentLength(URL url) {
        String contentLengthString = getHeaderField(url, "Content-Length");
        try {
            return Long.parseLong(contentLengthString);
        } catch (NumberFormatException e) {
            log.error("Error getting content length from: " + url.toString() + "\n" + "Content-length=" + contentLengthString);
            return 0;
        }

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

        if (url.getProtocol().toLowerCase().equals("ftp")) {
            return FTPUtils.resourceAvailable(url);
        }

        HttpURLConnection conn = null;
        try {
            conn = (HttpURLConnection) url.openConnection();
            int rc = conn.getResponseCode();
            return rc < 400;
        } catch (Exception e) {
            return false;
        }
        finally {
            if (conn != null) {
                conn.disconnect();
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

    public static class HttpClientInputStream extends FilterInputStream {
        GetMethod getMethod;

        protected HttpClientInputStream(GetMethod getMethod, InputStream inputStream) {
            super(inputStream);
            this.getMethod = getMethod;
        }

        @Override
        public void close() throws IOException {
            super.close();
            getMethod.releaseConnection();
        }
    }


}
