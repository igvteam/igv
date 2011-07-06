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
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.conn.params.ConnRoutePNames;
import org.apache.http.impl.client.BasicResponseHandler;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.io.*;
import java.net.URL;

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

    final static DefaultHttpClient client = new DefaultHttpClient();

    static {
        client.getParams().setParameter("http.protocol.allow-circular-redirects", true);
        client.getParams().setParameter("http.useragent", Globals.applicationString());
    }

    public static void setProxy(String proxyHost, int proxyPort, boolean auth, String user, String pw) {

        if (proxyHost != null) {
            HttpHost proxy = new HttpHost(proxyHost, proxyPort);
            client.getParams().setParameter(ConnRoutePNames.DEFAULT_PROXY, proxy);
            log.info("Proxy settings: " + proxyHost + ":" + proxyPort);
        }
        if (auth) {
            client.getCredentialsProvider().setCredentials(
                    new AuthScope(proxyHost, proxyPort),
                    new UsernamePasswordCredentials(user, pw));

        }
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
        HttpGet getMethod = new HttpGet(url.toExternalForm());
        HttpResponse response = execute(getMethod, url);
        return response.getEntity().getContent();
    }

    public static boolean resourceAvailable(URL url) {
        try {
            HttpHead headMethod = new HttpHead(url.toExternalForm());
            HttpResponse response = execute(headMethod, url);
            final int statusCode = response.getStatusLine().getStatusCode();
            return statusCode == 200;
        } catch (Exception e) {
            //log.error("Error checking for existence of resource: " + url);
            return false;
        }

    }

    /**
     * Execute a get on the url and return the header field value.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public static String getHeaderField(URL url, String key) throws IOException {
        HttpHead headMethod = new HttpHead(url.toExternalForm());
        HttpResponse response = execute(headMethod, url);
        return response.getFirstHeader(key).getValue();
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
    private static HttpResponse execute(HttpRequestBase getMethod, URL url) throws IOException {

        try {

            if (GSUtils.isGenomeSpace(url)) {
                GSUtils.checkForCookie(client, url);
            }
            HttpResponse response = client.execute(getMethod);
            final int statusCode = response.getStatusLine().getStatusCode();
            if (statusCode == 401) {
                // Try again                
                getMethod.abort();
                client.getCredentialsProvider().clear();
                login(url);
                return execute(getMethod, url);
            }
            if (statusCode != 200) {
                getMethod.abort();
                throw new RuntimeException("Error connecting.  Status code = " + statusCode);
            }
            return response;

        } catch (RuntimeException e) {
            // In case of an unexpected exception you may want to abort
            // the HTTP request in order to shut down the underlying
            // connection immediately. THis happens automatically for an IOException
            if (getMethod != null) getMethod.abort();
            throw e;
        }
    }


    private static void login(URL url) {

        Frame owner = IGV.hasInstance() ? IGV.getMainFrame() : null;

        // TODO -- only use the GS logo if this is a GS URL
        String userpass = getGSUserPass(owner);
        if (userpass == null) {
            throw new RuntimeException("Access denied:  " + url.toString());
        }
        UsernamePasswordCredentials GENOME_SPACE_CREDS = new UsernamePasswordCredentials(userpass);

        String host = GSUtils.isGenomeSpace(url) ? GSUtils.GENOME_SPACE_ID_SERVER : url.getHost();

        client.getCredentialsProvider().setCredentials(
                new AuthScope(host, AuthScope.ANY_PORT, AuthScope.ANY_REALM),
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


    public static String getGSUserPass(Frame owner) {

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

}
