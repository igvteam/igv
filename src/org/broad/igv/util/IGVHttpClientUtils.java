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

import org.apache.http.HttpEntity;
import org.apache.http.HttpHost;
import org.apache.http.HttpResponse;
import org.apache.http.auth.AuthScope;
import org.apache.http.auth.UsernamePasswordCredentials;
import org.apache.http.client.HttpClient;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.conn.params.ConnRoutePNames;
import org.apache.http.impl.client.BasicResponseHandler;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.http.impl.conn.tsccm.ThreadSafeClientConnManager;
import org.apache.http.util.EntityUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ftp.FTPClient;
import org.broad.igv.util.ftp.FTPStream;
import org.broad.igv.util.ftp.FTPUtils;
import org.broad.igv.util.stream.ApacheURLHelper;
import org.broad.tribble.util.SeekableHTTPStream;

import java.awt.*;
import java.io.*;
import java.net.InetSocketAddress;
import java.net.Proxy;
import java.net.URL;
import java.util.Map;

/**
 * New version of IGVHttpUtils built on Apache HttpClient 4.1.  Currently this version is only used for GenomeSpace
 * connections, which was its intention,  but eventually all client connections will use this class and IGVHttpUtils
 * will be eliminated.
 *
 * @author jrobinso
 * @date Jun 9, 2011
 */
public class IGVHttpClientUtils {

    private static Logger log = Logger.getLogger(IGVHttpClientUtils.class);

    final static DefaultHttpClient client;
    public static boolean byteRangeTested = false;
    public static boolean useByteRange = true;
    /**
     * Proxy settings (can be null)
     */
    //public static ProxySettings proxySettings = null;

    static {

        ThreadSafeClientConnManager cm = new ThreadSafeClientConnManager();
        cm.setMaxTotal(100);

        client = new DefaultHttpClient(cm);
        client.getParams().setParameter("http.protocol.allow-circular-redirects", true);
        client.getParams().setParameter("http.useragent", Globals.applicationString());
        
        
    }

    /**
     * Shutdown client and free all resources.  Called upon application exit.
     */
    public static void shutdown() {
        client.getConnectionManager().shutdown();
    }

    /**
     * Download the contents of the URL and save the results to a file.
     *
     * @throws IOException
     */
    public static boolean downloadFile(String url, File outputFile) throws IOException {

        log.info("Downloading " + url + " to " + outputFile.getAbsolutePath());

        HttpGet httpget = new HttpGet(url);

        HttpResponse response = client.execute(httpget);
        HttpEntity entity = response.getEntity();

        if (entity != null) {
            final long contentLength = entity.getContentLength();
            log.info("Content length = " + contentLength);

            InputStream is = null;
            OutputStream out = null;

            try {
                is = entity.getContent();
                out = new FileOutputStream(outputFile);

                byte[] buf = new byte[64 * 1024];
                int downloaded = 0;
                int bytesRead = 0;
                while ((bytesRead = is.read(buf)) != -1) {
                    out.write(buf, 0, bytesRead);
                    downloaded += bytesRead;
                }
                log.info("Download complete.  Total bytes downloaded = " + downloaded);
            }
            finally {
                if (is != null) is.close();
                if (out != null) {
                    out.flush();
                    out.close();
                }
            }
            long fileLength = outputFile.length();
            log.info("File length = " + fileLength);

            return contentLength <= 0 || contentLength == fileLength;


        }
        return false;

    }


    /**
     * Execute a get on the url and return the response stream.  It is the responsibility of the caller to
     * close the stream.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public static InputStream openConnectionStream(URL url) throws IOException {

        //TODO -- the protocol (ftp) test should be done before calling this method.
        if (url.getProtocol().toLowerCase().equals("ftp")) {
            String userInfo = url.getUserInfo();
            String host = url.getHost();
            String file = url.getPath();
            FTPClient ftp = FTPUtils.connect(host, userInfo);
            ftp.pasv();
            ftp.retr(file);
            return new FTPStream(ftp);
        } else {
            HttpGet getMethod = new HttpGet(url.toExternalForm());
            HttpResponse response = execute(getMethod, url);
            return new InputStreamWrapper(response);
        }
    }

    public static class InputStreamWrapper extends FilterInputStream {
        HttpResponse response;
        public InputStreamWrapper(HttpResponse response) throws IOException {
            super(response.getEntity().getContent());
            this.response = response;
        }

        @Override
        public void close() throws IOException {
            EntityUtils.consume(response.getEntity());
            super.close();
        }
    }


    /**
     * Test for the existence of a URL resource
     *
     * @param url
     * @return
     */
    public static boolean resourceAvailable(URL url) {
        try {
            HttpHead headMethod = new HttpHead(url.toExternalForm());
            HttpResponse response = execute(headMethod, url);
            final int statusCode = response.getStatusLine().getStatusCode();
            return statusCode == 200;
        } catch (Exception e) {
            //TODO -- non-expected exceptions (as opposed to "not found" exceptions) should be logged
            return false;
        }

    }

    /**
     * Execute a HEAD on the url and return the header field value.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public static String getHeaderField(URL url, String key) throws IOException {
        HttpHead headMethod = new HttpHead(url.toExternalForm());
        HttpResponse response = execute(headMethod, url);
        String value =  response.getFirstHeader(key).getValue();
        EntityUtils.consume(response.getEntity());
        return value;
    }

    /**
     * Execute a get request.  In the case of an authorization failure (401) this method is called recursively
     * with a login prompt until the correct credentials are entered,  or the user cancels triggering an
     * authorization exception.
     *
     * @param url
     * @return
     * @throws IOException
     */
    private static synchronized HttpResponse execute(HttpRequestBase method, URL url) throws IOException {

        try {

            if (GSUtils.isGenomeSpace(url)) {
                GSUtils.checkForCookie(client, url);
            }
            HttpResponse response = client.execute(method);
            final int statusCode = response.getStatusLine().getStatusCode();
            if (statusCode == 401) {
                // Try again                
                method.abort();
                client.getCredentialsProvider().clear();
                login(url);
                return execute(method, url);
            }
            if (statusCode != 200) {
                method.abort();
                throw new RuntimeException("Error connecting.  Status code = " + statusCode);
            }
            return response;

        } catch (RuntimeException e) {
            // An unexpected exception  -- abort the HTTP request in order to shut down the underlying
            // connection immediately. THis happens automatically for an IOException
            if (method != null) method.abort();
            throw e;
        }
    }

    /**
     * @param url
     */
    public static HttpResponse executeGet(URL url) {
        HttpGet get = new HttpGet(url.toExternalForm());
        try {
            return client.execute(get);
        } catch (IOException e) {
            log.error("Error executing get: " + url, e);
            return null;
        }

    }

    public static HttpResponse executeGet(URL url, Map<String, String> headers) {
        HttpGet get = new HttpGet(url.toExternalForm());
        for (Map.Entry<String, String> entry : headers.entrySet()) {
            get.setHeader(entry.getKey(), entry.getValue());
        }
        try {
            return client.execute(get);
        } catch (IOException e) {
            log.error("Error executing get: " + url, e);
            return null;
        }


    }

    private static void login(URL url) {

        Frame owner = IGV.hasInstance() ? IGV.getMainFrame() : null;

        String userpass = getUserPass(owner);
        if (userpass == null) {
            throw new RuntimeException("Access denied:  " + url.toString());
        }
        UsernamePasswordCredentials GENOME_SPACE_CREDS = new UsernamePasswordCredentials(userpass);

        String host = GSUtils.isGenomeSpace(url) ? GSUtils.GENOME_SPACE_ID_SERVER : url.getHost();

        client.getCredentialsProvider().setCredentials(
                new AuthScope(AuthScope.ANY_HOST, AuthScope.ANY_PORT, AuthScope.ANY_REALM),
                GENOME_SPACE_CREDS);

        if (GSUtils.isGenomeSpace(url)) {
            // Get the genomespace token
            try {
                HttpGet httpget = new HttpGet(GSUtils.identityServerUrl);
                ResponseHandler<String> responseHandler = new BasicResponseHandler();
                String responseBody = client.execute(httpget, responseHandler);
                if (responseBody != null && responseBody.length() > 0) {
                    String[] tokens = userpass.split(":");
                    String user = tokens[0];
                    GSUtils.saveLoginForSSO(responseBody, user);
                }
            } catch (IOException e) {
                log.error("Error fetching GS token", e);
            }

        }

    }

    /**
     * Open a modal login dialog and return
     *
     * @param owner
     * @return the user credentials in the form of "user:password".  If the  user cancels return null.
     */
    public static String getUserPass(Frame owner) {

        LoginDialog dlg = new LoginDialog(owner);
        dlg.setVisible(true);
        if (dlg.isCanceled()) {
            return null;
        } else {
            final String userString = dlg.getUsername();
            final String userPass = new String(dlg.getPassword());
            return userString + ":" + userPass;

        }

    }

    public static long getContentLength(URL url) {
        String contentLengthString = "";
        try {
            contentLengthString = getHeaderField(url, "Content-Length");
            return Long.parseLong(contentLengthString);
        } catch (IOException e) {
            log.error("Error getting content length from: " + url.toString(), e);
            return -1;

        }
        catch (NumberFormatException e) {
            log.error("Error getting content length from: " + url.toString() + "\n" + "Content-length=" + contentLengthString);
            return -1;
        }

    }

    public static long getContentLength(HttpResponse response) {
        String contentLengthString = "";
        try {
            contentLengthString = response.getFirstHeader("Content-Length").getValue();
            return Long.parseLong(contentLengthString);
        }
        catch (NumberFormatException e) {
            log.error("Error getting content length from: " + contentLengthString + "\n" + "Content-length=" + contentLengthString);
            return -1;
        }

    }


    public static boolean isURL(String string) {
        String lcString = string.toLowerCase();
        return lcString.startsWith("http://") || lcString.startsWith("https://") || lcString.startsWith("ftp://")
                || lcString.startsWith("file://");
    }

    public static boolean testByteRange() {

        try {
            String testURL = "http://www.broadinstitute.org/igvdata/byteRangeTest.txt";
            byte[] expectedBytes = {(byte) 'k', (byte) 'l', (byte) 'm', (byte) 'n', (byte) 'o'};

            SeekableHTTPStream str = new SeekableHTTPStream(new ApacheURLHelper(new URL(testURL)));
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

    public static boolean useByteRange() {
        useByteRange = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.USE_BYTE_RANGE);
        if (useByteRange && !byteRangeTested) {
            useByteRange = testByteRange();
            byteRangeTested = true;
        }
        return useByteRange;
    }


    public static void updateProxySettings() {

        String proxyHost;
        int proxyPort;
        boolean auth;
        String user;
        String pw = null;

        PreferenceManager prefMgr = PreferenceManager.getInstance();
        boolean useProxy = prefMgr.getAsBoolean(PreferenceManager.USE_PROXY);
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

        if (useProxy) {
            if (proxyHost != null) {
                HttpHost proxy = new HttpHost(proxyHost, proxyPort);
                client.getParams().setParameter(ConnRoutePNames.DEFAULT_PROXY, proxy);
                log.info("Proxy settings: " + proxyHost + ":" + proxyPort);
            }
            if (auth && pw != null) {
                client.getCredentialsProvider().setCredentials(
                        new AuthScope(proxyHost, proxyPort),
                        new UsernamePasswordCredentials(user, pw));
            }

        } else {
            client.getParams().removeParameter(ConnRoutePNames.DEFAULT_PROXY);
        }

    }

}
