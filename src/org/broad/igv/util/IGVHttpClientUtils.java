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
import org.apache.http.HttpHost;
import org.apache.http.HttpResponse;
import org.apache.http.auth.AuthScope;
import org.apache.http.auth.Credentials;
import org.apache.http.auth.UsernamePasswordCredentials;
import org.apache.http.auth.params.AuthPNames;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.client.methods.HttpPut;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.client.params.AuthPolicy;
import org.apache.http.conn.ClientConnectionManager;
import org.apache.http.conn.params.ConnRoutePNames;
import org.apache.http.conn.scheme.Scheme;
import org.apache.http.conn.scheme.SchemeRegistry;
import org.apache.http.conn.ssl.SSLSocketFactory;
import org.apache.http.entity.FileEntity;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.BasicResponseHandler;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.http.impl.conn.ProxySelectorRoutePlanner;
import org.apache.http.impl.conn.tsccm.ThreadSafeClientConnManager;
import org.apache.http.util.EntityUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;

import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ftp.FTPClient;
import org.broad.igv.util.ftp.FTPStream;
import org.broad.igv.util.ftp.FTPUtils;

import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
import java.awt.*;
import java.io.*;
import java.net.ProxySelector;
import java.net.URL;
import java.security.cert.CertificateException;
import java.security.cert.X509Certificate;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * New version of IGVHttpUtils built on Apache HttpClient 4.1.  Currently this version is only used for GenomeSpace
 * connections, which was its intention,  but eventually all client connections will use this class and IGVHttpUtils
 * will be eliminated.
 *
 * @author jrobinso
 * @date Jun 9, 2011
 */
public class IGVHttpClientUtils extends HttpUtils {

    private static Logger log = Logger.getLogger(IGVHttpClientUtils.class);

    private DefaultHttpClient client;
    private IdleConnectionMonitorThread monitorThread;

    private static IGVHttpClientUtils instance;

    static {
        synchronized (IGVHttpClientUtils.class) {
            instance = new IGVHttpClientUtils();
        }
    }

    public static IGVHttpClientUtils getInstance() {
        return instance;
    }

    private IGVHttpClientUtils() {
        client = createClient();
    }


    /**
     * Create the singleton client instance.  This is private to insure a single instance is created.
     *
     * @return
     */
    private DefaultHttpClient createClient() {

        try {
            ThreadSafeClientConnManager cm = new ThreadSafeClientConnManager();
            cm.setMaxTotal(100);

            SSLContext ctx = SSLContext.getInstance("TLS");
            X509TrustManager tm = getDefaultTrustManager();
            ctx.init(null, new TrustManager[]{tm}, null);
            SSLSocketFactory ssf = new SSLSocketFactory(ctx);
            ssf.setHostnameVerifier(SSLSocketFactory.ALLOW_ALL_HOSTNAME_VERIFIER);

            SchemeRegistry sr = cm.getSchemeRegistry();
            sr.register(new Scheme("https", ssf, 443));

            monitorThread = new IdleConnectionMonitorThread(cm);
            monitorThread.start();

            client = new DefaultHttpClient(cm);
            client.getParams().setParameter("http.protocol.allow-circular-redirects", true);
            client.getParams().setParameter("http.useragent", Globals.applicationString());


            boolean includeKerbeos = (new File(System.getenv("windir") + "\\krb5.ini").exists());
            ArrayList<String> authpref = new ArrayList<String>();
            authpref.add(AuthPolicy.BASIC);
            authpref.add(AuthPolicy.DIGEST);
            if (includeKerbeos) {
                authpref.add(AuthPolicy.SPNEGO);
            }
            authpref.add(AuthPolicy.NTLM);
            client.getParams().setParameter(AuthPNames.PROXY_AUTH_PREF, authpref);
            client.getParams().setParameter(AuthPNames.TARGET_AUTH_PREF, authpref);

            updateProxySettings();

            return client;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    /**
     * A default trust manager for SSL connections.  Basically trusts everybody.
     *
     * @return
     */
    private static X509TrustManager getDefaultTrustManager() {
        X509TrustManager tm = new X509TrustManager() {

            public void checkClientTrusted(X509Certificate[] xcs, String string) throws CertificateException {
            }

            public void checkServerTrusted(X509Certificate[] xcs, String string) throws CertificateException {
            }

            public X509Certificate[] getAcceptedIssuers() {
                return null;
            }
        };
        return tm;
    }

    /**
     * Shutdown client and free all resources.  Called upon application exit.
     */
    public void shutdown() {
        client.getConnectionManager().shutdown();
        monitorThread.shutdown();
    }

    /**
     * Return the contents of the url as a String.  This method should only be used for queries expected to return
     * a small amount of data.
     *
     * @param url
     * @return
     */
    public String getContentsAsString(URL url) throws IOException {
        HttpResponse response = executeGet(url);
        return EntityUtils.toString(response.getEntity());
    }


    /**
     * Execute a get on the url and return the response stream.  It is the responsibility of the caller to
     * close the stream.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public InputStream openConnectionStream(URL url) throws IOException {
        return openConnectionStream(url, false);
    }


    /**
     * Execute a get on the url and return the response stream.  Optionally abort httpget when closing the stream.
     * This should be done when reading only a portion of the response stream, as otherwise the close() method in the
     * HttpClient stream class will read the rest of the stream before closing it.
     * <p/>
     * <p/>
     * NOTE:  A better solution for the partial read problem would be to use byte-range headers to get only the portion
     * of the file needed.
     * <p/>
     * It is the responsibility of the caller to
     * close the stream.
     *
     * @param url
     * @param abortOnClose true if HttpGet.abort() should be called upon close.  Note this will also kill the connection.
     * @return
     * @throws IOException
     */

    public InputStream openConnectionStream(URL url, boolean abortOnClose) throws IOException {

        return openConnectionStream(url, abortOnClose, null);
    }


    public InputStream openConnectionStream(URL url, Map<String, String> headers) throws IOException {
        HttpGet get = new HttpGet(url.toExternalForm());
        for (Map.Entry<String, String> entry : headers.entrySet()) {
            get.setHeader(entry.getKey(), entry.getValue());
        }
        return execute(get, url).getEntity().getContent();
    }


    public InputStream openConnectionStream(URL url, boolean abortOnClose, Map<String, String> headers) throws IOException {

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
            if (headers != null) {
                for (Map.Entry<String, String> entry : headers.entrySet()) {
                    getMethod.setHeader(entry.getKey(), entry.getValue());
                }
            }

            HttpResponse response = execute(getMethod, url);

            // Wrap the response stream to do extra cleanaup upon close.
            InputStream is = response.getEntity().getContent();
            return abortOnClose ? new ResponseInputStream(getMethod, is) : is;

        }
    }


    /**
     * Test for the existence of a URL resource
     *
     * @param url
     * @return
     */
    public boolean resourceAvailable(URL url) {

        try {
            if (GSUtils.isGenomeSpace(url.toExternalForm())) {
                //The GenomeSpace server does not support "HEAD".  DO a get and abort the method
                //after retrieving the status code
                HttpGet getMethod = new HttpGet(url.toExternalForm());
                HttpResponse response = execute(getMethod, url);
                int statusCode = response.getStatusLine().getStatusCode();
                getMethod.abort();
                return statusCode == 200;

            } else {
                HttpHead headMethod = new HttpHead(url.toExternalForm());
                HttpResponse response = execute(headMethod, url);
                final int statusCode = response.getStatusLine().getStatusCode();
                return statusCode == 200;
            }
        } catch (FileNotFoundException e) {
            return false;
        } catch (Exception e) {
            log.error("Error checking resoruce availability", e);
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
    public String getHeaderField(URL url, String key) throws IOException {

        if (GSUtils.isGenomeSpace(url.toExternalForm())) {
            //The GenomeSpace server does not support "HEAD".  DO a get and abort the method
            //after retrieving the header value
            HttpGet getMethod = new HttpGet(url.toExternalForm());
            HttpResponse response = execute(getMethod, url);
            String value = null;
            Header header = response.getFirstHeader(key);
            if (header != null) {
                value = header.getValue();
            }
            getMethod.abort();
            return value;

        } else {
            HttpHead headMethod = new HttpHead(url.toExternalForm());
            HttpResponse response = execute(headMethod, url);

            String value = null;
            Header header = response.getFirstHeader(key);
            if (header != null) {
                value = header.getValue();
            }
            return value;
        }
    }

    /**
     * @param url
     */
    private HttpResponse executeGet(URL url) throws IOException {
        HttpGet get = new HttpGet(url.toExternalForm());
        return execute(get, url);

    }


    /**
     * Upload a file.
     * <p/>
     * Note: this method was written for, and has only been tested against, the GenomeSpace amazon server.
     *
     * @throws IOException
     */
    public String createGenomeSpaceDirectory(URL url, String body) throws IOException {
        HttpPut put = new HttpPut(url.toExternalForm());
        put.setHeader("Content-Type", "application/json");
        StringEntity se = new StringEntity(body);
        put.setEntity(se);

        HttpResponse response = execute(put, url);
        String responseString = EntityUtils.toString(response.getEntity());

        int code = response.getStatusLine().getStatusCode();
        if (code != 200) {
            throw new IOException("Error creating directory: " + code);
        }

        return responseString;


    }

    /**
     * Upload a file.
     * <p/>
     * Note: this method was written for, and has only been tested against, the GenomeSpace amazon server.
     *
     *
     *
     * @param uri
     * @param file
     * @param headers
     * @throws IOException
     */
    public void uploadGenomeSpaceFile(String uri, File file, Map<String, String> headers) throws IOException {

        HttpPut put = new HttpPut(uri);
        try {
            FileEntity entity = new FileEntity(file, "text");

            put.setEntity(entity);
            if (headers != null) {
                for (Map.Entry<String, String> entry : headers.entrySet()) {
                    put.addHeader(entry.getKey(), entry.getValue());
                }
            }

            HttpResponse response = client.execute(put);
            EntityUtils.consume(response.getEntity());

            final int statusCode = response.getStatusLine().getStatusCode();
            if (statusCode == 401) {
                // Try again
                client.getCredentialsProvider().clear();
                login(new URL(uri));
                uploadGenomeSpaceFile(uri, file, headers);
            } else if (statusCode == 404 || statusCode == 410) {
                put.abort();
                throw new FileNotFoundException();
            } else if (statusCode >= 400) {
                put.abort();
                throw new HttpResponseException(statusCode);
            }

        } catch (RuntimeException e) {
            // An unexpected exception  -- abort the HTTP request in order to shut down the underlying
            // connection immediately. THis happens automatically for an IOException
            if (put != null) put.abort();
            throw e;
        }

    }

    /**
     * Execute a  request.  In the case of an authorization failure (401) this method is called recursively
     * with a login prompt until the correct credentials are entered,  or the user cancels triggering an
     * authorization exception.
     *
     * @param url
     * @return
     * @throws IOException
     */
    private HttpResponse execute(HttpRequestBase method, URL url) throws IOException {

        try {

            if (GSUtils.isGenomeSpace(url.toString())) {
                GSUtils.checkForCookie(client, url.getHost());
            }

            HttpResponse response = client.execute(method);
            final int statusCode = response.getStatusLine().getStatusCode();
            if (statusCode == 401) {
                // Try again
                client.getCredentialsProvider().clear();
                login(url);
                return execute(method, url);
            } else if (statusCode == 404 || statusCode == 410) {
                method.abort();
                throw new FileNotFoundException();
            } else if (statusCode >= 400) {
                method.abort();
                throw new HttpResponseException(statusCode);
            }
            return response;

        } catch (RuntimeException e) {
            // An unexpected exception  -- abort the HTTP request in order to shut down the underlying
            // connection immediately. THis happens automatically for an IOException
            if (method != null) method.abort();
            throw e;
        }
    }


    private void login(URL url) {

        String server = url.getHost();

        Frame owner = IGV.hasInstance() ? IGV.getMainFrame() : null;
        final boolean genomeSpace = GSUtils.isGenomeSpace(server);

        String userpass = getUserPass(owner, url, genomeSpace);
        if (userpass == null) {
            throw new RuntimeException("Access denied:  " + server);
        }
        UsernamePasswordCredentials GENOME_SPACE_CREDS = new UsernamePasswordCredentials(userpass);

        client.getCredentialsProvider().setCredentials(
                new AuthScope(AuthScope.ANY_HOST, AuthScope.ANY_PORT, AuthScope.ANY_REALM),
                GENOME_SPACE_CREDS);

        if (genomeSpace) {
            // Get the genomespace token
            HttpGet httpget = null;
            try {
                httpget = new HttpGet(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ID_SERVER));
                ResponseHandler<String> responseHandler = new BasicResponseHandler();
                String responseBody = client.execute(httpget, responseHandler);
                if (responseBody != null && responseBody.length() > 0) {
                    String[] tokens = userpass.split(":");
                    String user = tokens[0];
                    GSUtils.saveGSLogin(responseBody, user);
                }
            } catch (org.apache.http.client.HttpResponseException e) {
                int code = e.getStatusCode();
                if(code == 401) {
                    MessageUtils.showMessage("Invalid username or password.");
                }
                if (httpget != null) {
                    httpget.abort();
                }
                createClient();
            } catch (IOException e) {
                log.error("Error fetching GS token", e);
                if (httpget != null) {
                    httpget.abort();
                }
                createClient();

            }

        }

    }

    /**
     * Remove the GenomeSpace cookie
     */
    public void removeGSCookie() {
        createClient();
//        CookieStore cookieStore = client.getCookieStore();
//        // Copy the list
//        java.util.List<Cookie> cookies = cookieStore.getCookies();
//        cookieStore.clear();
//        for (Cookie cookie : cookies) {
//            if (!cookie.getName().equals(GSUtils.AUTH_TOKEN_COOKIE_NAME)) {
//                cookieStore.addCookie(cookie);
//            }
//        }
//
//        client.getCredentialsProvider().clear();

    }


    /**
     * Open a modal login dialog and return
     *
     * @param owner
     * @param url
     * @return the user credentials in the form of "user:password".  If the  user cancels return null.
     */
    public static String getUserPass(Frame owner, URL url, boolean isGenomeSpace) {


        LoginDialog dlg = new LoginDialog(owner, isGenomeSpace, url.toString(), false);
        //  GSUtils.isGenomeSpace(server)
        // public LoginDialog(Frame owner, boolean isGenomeSpace, String resource, boolean proxyChallenge) {
        dlg.setVisible(true);
        if (dlg.isCanceled()) {
            return null;
        } else {
            final String userString = dlg.getUsername();
            final String userPass = new String(dlg.getPassword());
            return userString + ":" + userPass;

        }

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
                Credentials creds = new UsernamePasswordCredentials(user, pw);
                client.getCredentialsProvider().setCredentials(new AuthScope(proxyHost, proxyPort), creds);
            }

        } else {
            client.getParams().removeParameter(ConnRoutePNames.DEFAULT_PROXY);
            ProxySelectorRoutePlanner routePlanner = new ProxySelectorRoutePlanner(
                    client.getConnectionManager().getSchemeRegistry(),
                    ProxySelector.getDefault());
            client.setRoutePlanner(routePlanner);
        }

    }

    /**
     * Download the contents of the URL and save the results to a file.
     *
     * @throws IOException
     */
    public boolean downloadFile(String url, File outputFile) throws IOException {

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
            } catch (IOException e) {
                httpget.abort();
                throw e;
            } finally {
                if (is != null) is.close();
                if (out != null) {
                    out.flush();
                    out.close();
                }
            }
            long fileLength = outputFile.length();
            //log.info("File length = " + fileLength);

            return contentLength <= 0 || contentLength == fileLength;


        }
        return false;

    }


    /**
     * Wrapper for an Entity input stream.
     * <p/>
     * NOTE: Without the call to getMethod.abort() in the close method the entire contents of the underlying
     * stream will be read by the Apache implementation.  As IGV peeks at large files by opening a stream and
     * reading a few lines we have to protect against this, or the application might hang.  For the future -- a better
     * solution would be to use byte-range requests to request only digestible sections of the file at a time.
     */
    public static class ResponseInputStream extends FilterInputStream {
        HttpGet getMethod;


        public ResponseInputStream(HttpGet getMethod, InputStream content) {
            super(content);
            this.getMethod = getMethod;
        }

        @Override
        public void close() throws IOException {
            getMethod.abort();
        }
    }


    /**
     * Thread to flush idle connections periodically
     */
    public static class IdleConnectionMonitorThread extends Thread {

        private final ClientConnectionManager connMgr;
        private volatile boolean shutdown;

        public IdleConnectionMonitorThread(ClientConnectionManager connMgr) {
            super();
            this.connMgr = connMgr;
        }

        @Override
        public void run() {
            try {
                while (!shutdown) {
                    synchronized (this) {
                        wait(60000);
                        // Close expired connections
                        connMgr.closeExpiredConnections();
                        // Optionally, close connections
                        // that have been idle longer than 300 sec
                        connMgr.closeIdleConnections(300, TimeUnit.SECONDS);
                    }
                }
            } catch (InterruptedException ex) {
                // terminate
            }
        }

        public void shutdown() {
            shutdown = true;
            synchronized (this) {
                notifyAll();
            }
        }

    }


    /**
     * Returns the proxy information for the specified sampleURL using JRE 1.4
     * specific plugin classes.
     *
     * Notes:
     *     Plugin 1.4 Final added
     *     com.sun.java.browser.net.* classes ProxyInfo & ProxyService...
     *     Use those with JREs => 1.4
     *
     * @param sampleURL the URL to check proxy settings for
     * @return ProxyHost the host and port of the proxy that should be used
     */
    /* private static ProxyHost detectProxySettingsJDK14_JDK15_JDK16(URL sampleURL) {
      ProxyHost result = null;
      try {
          // Look around for the 1.4.X plugin proxy detection class...
          // Without it, cannot autodetect...
          Class ProxyServiceClass =
              Class.forName("com.sun.java.browser.net.ProxyService");
          Method getProxyInfoMethod =
              ProxyServiceClass.getDeclaredMethod("getProxyInfo",
                                                  new Class[] {URL.class});
          Object proxyInfoArrayObj =
              getProxyInfoMethod.invoke(null, new Object[] {sampleURL});

          if (proxyInfoArrayObj == null
                  || Array.getLength(proxyInfoArrayObj) == 0) {
              if (log.isDebugEnabled()) {
                  log.debug("1.4.X reported NULL proxy (no proxy assumed)");
              }
              result = NO_PROXY_HOST;
          } else {
              Object proxyInfoObject = Array.get(proxyInfoArrayObj, 0);
              Class proxyInfoClass = proxyInfoObject.getClass();
              Method getHostMethod =
                  proxyInfoClass.getDeclaredMethod("getHost",null);
              String proxyIP =
                  (String)getHostMethod.invoke(proxyInfoObject, null);
              Method getPortMethod =
                  proxyInfoClass.getDeclaredMethod("getPort",null);
              Integer portInteger =
                  (Integer)getPortMethod.invoke(proxyInfoObject, null);
              int proxyPort = portInteger.intValue();
              if (log.isDebugEnabled()) {
                  log.debug("1.4.X Proxy info geProxy:"+proxyIP+
                            " get Port:"+proxyPort);
              }
              result = new ProxyHost(proxyIP, proxyPort);
          }
      } catch (Exception e) {
          e.printStackTrace();
          log.warn("Sun Plugin 1.4.X proxy detection class not found, " +
                   "will try failover detection, e:"+e);
      }
      return result;
  }  */

}
