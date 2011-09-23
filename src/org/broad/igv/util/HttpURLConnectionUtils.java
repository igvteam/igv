/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

import org.apache.http.Header;
import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpHead;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ftp.FTPClient;
import org.broad.igv.util.ftp.FTPStream;
import org.broad.igv.util.ftp.FTPUtils;
import sun.misc.BASE64Encoder;

import java.awt.*;
import java.io.*;
import java.net.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Notes -- 401 => client authentication,  407 => proxy authentication,  403 => forbidden
 *
 * @author Jim Robinson
 * @date 9/22/11
 */
public class HttpURLConnectionUtils extends HttpUtils {

    private static Logger log = Logger.getLogger(IGVHttpClientUtils.class);

    private static ProxySettings proxySettings = null;

    private static HttpURLConnectionUtils instance;

    static {
        synchronized (HttpURLConnectionUtils.class) {
            instance = new HttpURLConnectionUtils();
        }
    }

    public static HttpURLConnectionUtils getInstance() {
        return instance;
    }

    private HttpURLConnectionUtils() {
        Authenticator.setDefault(new IGVAuthenticator());
    }

    public void shutdown() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
      * Return the contents of the url as a String.  This method should only be used for queries expected to return
      * a small amount of data.
      *
      * @param url
      * @return
      */
    public String getContentsAsString(URL url) throws IOException {
        return openConnection(url, null).getContent().toString();
    }

    public InputStream openConnectionStream(URL url) throws IOException {
        if (url.getProtocol().toUpperCase().equals("FTP")) {
            String userInfo = url.getUserInfo();
            String host = url.getHost();
            String file = url.getPath();
            FTPClient ftp = FTPUtils.connect(host, userInfo);
            ftp.pasv();
            ftp.retr(file);
            return new FTPStream(ftp);

        } else {
            return openConnectionStream(url, null);
        }
    }

    public InputStream openConnectionStream(URL url, boolean abortOnClose) throws IOException {
        return openConnection(url, null).getInputStream();
    }

    public InputStream openConnectionStream(URL url, Map<String, String> requestProperties) throws IOException {

        HttpURLConnection conn = openConnection(url, requestProperties);
        return conn.getInputStream();

    }

    public InputStream openConnectionStream(URL url, boolean abortOnClose, Map<String, String> requestProperties) throws IOException {
        HttpURLConnection conn = openConnection(url, requestProperties);
        return conn.getInputStream();
    }

    public boolean resourceAvailable(URL url) {

        try {
            HttpURLConnection conn = openConnection(url, null, "HEAD");
            int code = conn.getResponseCode();
            return code == 200;
        } catch (IOException e) {
            return false;
        }
    }

    public String getHeaderField(URL url, String key) throws IOException {
        HttpURLConnection conn = openConnection(url, null, "HEAD");
        int code = conn.getResponseCode();
        // TODO -- check code
        return conn.getHeaderField(key);
    }

    public long getContentLength(URL url) throws IOException {
        String contentLengthString = "";

        contentLengthString = getHeaderField(url, "Content-Length");
        if (contentLengthString == null) {
            return -1;
        } else {
            return Long.parseLong(contentLengthString);
        }
    }

    public void updateProxySettings() {
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


         proxySettings = new ProxySettings(useProxy, user, pw, auth, proxyHost, proxyPort);
    }

    public boolean downloadFile(String url, File outputFile) throws IOException {

        log.info("Downloading " + url + " to " + outputFile.getAbsolutePath());

        HttpURLConnection conn = openConnection(new URL(url), null);
        int code = conn.getResponseCode();
        // TODO -- check code

        long contentLength = -1;
        String contentLengthString = conn.getHeaderField("Content-Length");
        if (contentLengthString != null) {
            contentLength = Long.parseLong(contentLengthString);
        }


        log.info("Content length = " + contentLength);

        InputStream is = null;
        OutputStream out = null;

        try {
            is = conn.getInputStream();
            out = new FileOutputStream(outputFile);

            byte[] buf = new byte[64 * 1024];
            int downloaded = 0;
            int bytesRead = 0;
            while ((bytesRead = is.read(buf)) != -1) {
                out.write(buf, 0, bytesRead);
                downloaded += bytesRead;
            }
            log.info("Download complete.  Total bytes downloaded = " + downloaded);
        } finally {
            if (is != null) is.close();
            if (out != null) {
                out.flush();
                out.close();
            }
        }
        long fileLength = outputFile.length();

        return contentLength <= 0 || contentLength == fileLength;
    }

    @Override
    public void uploadFile(URI uri, File localFile, Map<String, String> headers) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String createDirectory(URL url, String body) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    private static HttpURLConnection openConnection(URL url, Map<String, String> requestProperties) throws IOException {
        return openConnection(url, requestProperties, "GET");
    }

    private static HttpURLConnection openConnection(URL url, Map<String, String> requestProperties, String method) throws IOException {

        boolean useProxy = proxySettings != null && proxySettings.useProxy && proxySettings.proxyHost != null &&
                proxySettings.proxyPort > 0;

        HttpURLConnection conn = null;
        if (useProxy) {
            Proxy proxy = new Proxy(Proxy.Type.HTTP, new InetSocketAddress(proxySettings.proxyHost, proxySettings.proxyPort));
            conn = (HttpURLConnection) url.openConnection(proxy);

            if (proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                byte[] bytes = (proxySettings.user + ":" + proxySettings.pw).getBytes();
                String encodedUserPwd = (new BASE64Encoder()).encode(bytes);
                conn.setRequestProperty("Proxy-Authorization", "Basic " + encodedUserPwd);
            }
        } else {
            conn = (HttpURLConnection) url.openConnection();

        }
        conn.setConnectTimeout(10000);
        conn.setReadTimeout(60000);
        conn.setRequestMethod(method);
        conn.setRequestProperty("Connection", "close");
        if (requestProperties != null) {
            for (Map.Entry<String, String> prop : requestProperties.entrySet()) {
                conn.setRequestProperty(prop.getKey(), prop.getValue());
            }
        }
        return conn;
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

    public static class IGVAuthenticator extends Authenticator {

        @Override
        protected PasswordAuthentication getPasswordAuthentication() {

            Authenticator.RequestorType type = getRequestorType();
            if (type == RequestorType.PROXY) {
                // Todo -- use stored proxy user/pw?
            }

            Frame owner = IGV.hasInstance() ? IGV.getMainFrame() : null;
            LoginDialog dlg = new LoginDialog(owner);
            dlg.setVisible(true);
            if (dlg.isCanceled()) {
                return null;
            } else {
                final String userString = dlg.getUsername();
                final char[] userPass = dlg.getPassword();
                return new PasswordAuthentication(userString, userPass);

            }
        }
    }

}
