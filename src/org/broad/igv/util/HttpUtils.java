/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
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

import biz.source_code.base64Coder.Base64Coder;
import org.apache.log4j.Logger;
import org.apache.tomcat.util.HttpDate;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.stream.IGVUrlHelper;
import org.broad.tribble.util.SeekableHTTPStream;
import org.broad.tribble.util.ftp.FTPClient;
import org.broad.tribble.util.ftp.FTPStream;
import org.broad.tribble.util.ftp.FTPUtils;

import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
import java.awt.*;
import java.io.*;
import java.net.*;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Wrapper utility class... for interacting with HttpURLConnection.
 *
 * @author Jim Robinson
 * @date 9/22/11
 */
public class HttpUtils {

    private static Logger log = Logger.getLogger(HttpUtils.class);

    public static boolean byteRangeTested = false;
    public static boolean byteRangeTestSuccess = false;
    private static HttpUtils instance;

    private ProxySettings proxySettings = null;
    private final int MAX_REDIRECTS = 5;


    /**
     *  Create the single instance  and register the cookie manager
     */
    static {
        synchronized (HttpUtils.class) {
            org.broad.tribble.util.ParsingUtils.registerHelperClass(IGVUrlHelper.class);
            instance = new HttpUtils();
            instance.disableCertificateValidation();
            CookieHandler.setDefault(new CookieManager());
        }
    }

    /**
     * @return the single instance
     */
    public static HttpUtils getInstance() {
        return instance;
    }

    /**
     * Constructor
     */
    private HttpUtils() {
        Authenticator.setDefault(new IGVAuthenticator());
    }

    public static boolean isURL(String string) {
        String lcString = string.toLowerCase();
        return lcString.startsWith("http://") || lcString.startsWith("https://") || lcString.startsWith("ftp://")
                || lcString.startsWith("file://");
    }

    /**
     * Join the {@code elements} with the character {@code joiner},
     * URLencoding the {@code elements} along the way. {@code joiner}
     * is NOT URLEncoded
     * Example:
     * String[] parm_list = new String[]{"app les", "oranges", "bananas"};
     * String formatted = buildURLString(Arrays.asList(parm_list), "+");
     * <p/>
     * formatted will be "app%20les+oranges+bananas"
     *
     * @param elements
     * @param joiner
     * @return
     */
    public static String buildURLString(Iterable<String> elements, String joiner) {

        Iterator<String> iter = elements.iterator();
        if (!iter.hasNext()) {
            return "";
        }
        String wholequery = iter.next();
        try {
            while (iter.hasNext()) {
                wholequery += joiner + URLEncoder.encode(iter.next(), "UTF-8");
            }
            return wholequery;
        } catch (UnsupportedEncodingException e) {
            throw new IllegalArgumentException("Bad argument in genelist: " + e.getMessage());
        }
    }

    /**
     * Test to see if this client can successfully retrieve a portion of a remote file using the byte-range header.
     * This is not a test of the server, but the client.  In some environments the byte-range header gets removed
     * by filters after the request is made by IGV.
     *
     * @return
     */
    public static boolean testByteRange() {

        log.info("Testing range-byte request");
        try {
            String testURL = "http://www.broadinstitute.org/igvdata/annotations/seq/hg19/chr12.txt";
            byte[] expectedBytes = {'T', 'C', 'G', 'C', 'T', 'T', 'G', 'A', 'A', 'C', 'C', 'C', 'G', 'G',
                    'G', 'A', 'G', 'A', 'G', 'G'};


            SeekableHTTPStream str = new SeekableHTTPStream(new IGVUrlHelper(new URL(testURL)));
            str.seek(25350000);
            byte[] buffer = new byte[80000];
            str.read(buffer);

//            for(int i=0; i<expectedBytes.length; i++) {
//
//            }

            for (int i = 0; i < expectedBytes.length; i++) {
                if (buffer[i] != expectedBytes[i]) {
                    log.info("Range-byte test failed -- problem with client network environment.");
                    return false;
                }
            }
            log.info("Range-byte request succeeded");
            return true;
        } catch (IOException e) {
            log.error("Error while testing byte range " + e.getMessage());
            // We could not reach the test server, so we can't know if this client can do byte-range tests or
            // not.  Take the "optimistic" view.
            return true;
        }
    }



    public static boolean useByteRange(URL url) {


        // We can test byte-range success for broad hosted data. We can't know if they work or not in other
        // environments (e.g. intranets)
        if (url.getHost().contains("broadinstitute.org")) {
            // Test broad urls for successful byte range requests.
            if (!byteRangeTested) {
                byteRangeTestSuccess = testByteRange();
                byteRangeTested = true;   // <= to prevent testing again

            }
            return byteRangeTestSuccess;
        } else {
            return true;
        }


    }

    public void shutdown() {
        // Do any cleanup required here
    }


    /**
     * Return the contents of the url as a String.  This method should only be used for queries expected to return
     * a small amount of data.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public String getContentsAsString(URL url) throws IOException {

        InputStream is = null;
        HttpURLConnection conn = openConnection(url, null);
        try {
            is = conn.getInputStream();
            return readContents(is);

        } finally {
            if (is != null) is.close();
        }
    }

    public String getContentsAsJSON(URL url) throws IOException {

        InputStream is = null;
        Map<String, String> reqProperties = new HashMap();
        reqProperties.put("Accept", "application/json,text/plain");
        HttpURLConnection conn = openConnection(url, reqProperties);
        try {
            is = conn.getInputStream();
            return readContents(is);

        } finally {
            if (is != null) is.close();
        }
    }

    /**
     * Open a connection stream for the URL.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public InputStream openConnectionStream(URL url) throws IOException {
        if (url.getProtocol().toLowerCase().equals("ftp")) {
            String userInfo = url.getUserInfo();
            String host = url.getHost();
            String file = url.getPath();
            FTPClient ftp = FTPUtils.connect(host, userInfo, new UserPasswordInputImpl());
            ftp.pasv();
            ftp.retr(file);
            return new FTPStream(ftp);

        } else {
            return openConnectionStream(url, null);
        }
    }

    public InputStream openConnectionStream(URL url, Map<String, String> requestProperties) throws IOException {

        HttpURLConnection conn = openConnection(url, requestProperties);
        InputStream input = conn.getInputStream();
        if ("gzip".equals(conn.getContentEncoding())) {
            input = new GZIPInputStream(input);
        }
        return input;
    }


    public boolean resourceAvailable(URL url) {

        if (url.getProtocol().toLowerCase().equals("ftp")) {
            return FTPUtils.resourceAvailable(url);
        }

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
        return conn.getHeaderField(key);
    }

    public long getContentLength(URL url) throws IOException {

        String contentLengthString = getHeaderField(url, "Content-Length");
        if (contentLengthString == null) {
            return -1;
        } else {
            return Long.parseLong(contentLengthString);
        }
    }

    /**
     * Compare a local and remote resource.
     *
     * @param file
     * @param url
     * @return true if the files are "the same", false if the remote file has been modified wrt the local one.
     * @throws IOException
     */
    public boolean compareResources(File file, URL url) throws IOException {


        if (!file.exists()) {
            return false;
        }

        HttpURLConnection conn = openConnection(url, null, "HEAD");

        // Check content-length first
        long contentLength = -1;
        String contentLengthString = conn.getHeaderField("Content-Length");
        if (contentLengthString != null) {
            try {
                contentLength = Long.parseLong(contentLengthString);
            } catch (NumberFormatException e) {
                log.error("Error parsing content-length string: " + contentLengthString + " from URL: "
                        + url.toString());
                contentLength = -1;
            }
        }
        if (contentLength != file.length()) {
            return false;
        }

        // Compare last-modified dates
        String lastModifiedString = conn.getHeaderField("Last-Modified");
        if (lastModifiedString == null) {
            return false;
        } else {
            HttpDate date = new HttpDate();
            date.parse(lastModifiedString);
            long remoteModifiedTime = date.getTime();
            long localModifiedTime = file.lastModified();
            return remoteModifiedTime <= localModifiedTime;
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


    public void uploadGenomeSpaceFile(String uri, File file, Map<String, String> headers) throws IOException {

        HttpURLConnection urlconnection = null;
        OutputStream bos = null;

        URL url = new URL(uri);
        urlconnection = openConnection(url, headers, "PUT");
        urlconnection.setDoOutput(true);
        urlconnection.setDoInput(true);

        bos = new BufferedOutputStream(urlconnection.getOutputStream());
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(file));
        int i;
        // read byte by byte until end of stream
        while ((i = bis.read()) > 0) {
            bos.write(i);
        }
        bos.close();
        int responseCode = urlconnection.getResponseCode();

        // Error messages below.
        if (responseCode >= 400) {
            String message = readErrorStream(urlconnection);
            throw new IOException("Error uploading " + file.getName() + " : " + message);
        }
    }


    public String createGenomeSpaceDirectory(URL url, String body) throws IOException {

        HttpURLConnection urlconnection = null;
        OutputStream bos = null;

        Map<String, String> headers = new HashMap<String, String>();
        headers.put("Content-Type", "application/json");
        headers.put("Content-Length", String.valueOf(body.getBytes().length));

        urlconnection = openConnection(url, headers, "PUT");
        urlconnection.setDoOutput(true);
        urlconnection.setDoInput(true);

        bos = new BufferedOutputStream(urlconnection.getOutputStream());
        bos.write(body.getBytes());
        bos.close();
        int responseCode = urlconnection.getResponseCode();

        // Error messages below.
        StringBuffer buf = new StringBuffer();
        InputStream inputStream;

        if (responseCode >= 200 && responseCode < 300) {
            inputStream = urlconnection.getInputStream();
        } else {
            inputStream = urlconnection.getErrorStream();
        }
        BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
        String nextLine;
        while ((nextLine = br.readLine()) != null) {
            buf.append(nextLine);
            buf.append('\n');
        }
        inputStream.close();

        if (responseCode >= 200 && responseCode < 300) {
            return buf.toString();
        } else {
            throw new IOException("Error creating GS directory: " + buf.toString());
        }
    }

    /**
     * Code for disabling SSL certification
     */
    private void disableCertificateValidation() {
        // Create a trust manager that does not validate certificate chains
        TrustManager[] trustAllCerts = new TrustManager[]{
                new X509TrustManager() {
                    public java.security.cert.X509Certificate[] getAcceptedIssuers() {
                        return null;
                    }

                    public void checkClientTrusted(
                            java.security.cert.X509Certificate[] certs, String authType) {
                    }

                    public void checkServerTrusted(
                            java.security.cert.X509Certificate[] certs, String authType) {
                    }
                }
        };

        // Install the all-trusting trust manager
        try {
            SSLContext sc = SSLContext.getInstance("SSL");
            sc.init(null, trustAllCerts, null);
            HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
        } catch (NoSuchAlgorithmException e) {
        } catch (KeyManagementException e) {
        }

    }

    private String readContents(InputStream is) throws IOException {
        BufferedInputStream bis = new BufferedInputStream(is);
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        int b;
        while ((b = bis.read()) >= 0) {
            bos.write(b);
        }
        return new String(bos.toByteArray());
    }

    private String readErrorStream(HttpURLConnection connection) throws IOException {
        InputStream inputStream = null;

        try {
            inputStream = connection.getErrorStream();
            if (inputStream == null) {
                return null;
            }
            return readContents(inputStream);
        } finally {
            if (inputStream != null) inputStream.close();
        }


    }

    private HttpURLConnection openConnection(URL url, Map<String, String> requestProperties) throws IOException {
        return openConnection(url, requestProperties, "GET");
    }

    private HttpURLConnection openConnection(URL url, Map<String, String> requestProperties, String method) throws IOException {
        return openConnection(url, requestProperties, method, 0);
    }

    /**
     * The "real" connection method
     *
     * @param url
     * @param requestProperties
     * @param method
     * @return
     * @throws java.io.IOException
     */
    private HttpURLConnection openConnection(
            URL url, Map<String, String> requestProperties, String method, int redirectCount) throws IOException {

        boolean useProxy = proxySettings != null && proxySettings.useProxy && proxySettings.proxyHost != null &&
                proxySettings.proxyPort > 0;

        HttpURLConnection conn;
        if (useProxy) {
            Proxy proxy = new Proxy(Proxy.Type.HTTP, new InetSocketAddress(proxySettings.proxyHost, proxySettings.proxyPort));
            conn = (HttpURLConnection) url.openConnection(proxy);

            if (proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                byte[] bytes = (proxySettings.user + ":" + proxySettings.pw).getBytes();

                String encodedUserPwd = String.valueOf(Base64Coder.encode(bytes));
                conn.setRequestProperty("Proxy-Authorization", "Basic " + encodedUserPwd);
            }
        } else {
            conn = (HttpURLConnection) url.openConnection();
        }

        if (GSUtils.isGenomeSpace(url)) {

            String token = GSUtils.getGSToken();
            if (token != null) conn.setRequestProperty("Cookie", "gs-token=" + token);
            //conn.setRequestProperty("Accept", "application/json,text/plain");

        }

        conn.setConnectTimeout(Globals.CONNECT_TIMEOUT);
        conn.setReadTimeout(Globals.READ_TIMEOUT);
        conn.setRequestMethod(method);
        conn.setRequestProperty("Connection", "Keep-Alive");
        if (requestProperties != null) {
            for (Map.Entry<String, String> prop : requestProperties.entrySet()) {
                conn.setRequestProperty(prop.getKey(), prop.getValue());
            }
        }
        conn.setRequestProperty("User-Agent", Globals.applicationString());

        if (method.equals("PUT")) {
            return conn;
        } else {
            int code = conn.getResponseCode();
            if (code >= 200 && code < 300) {
                // A genome-space hack.  We want to catch redirects from the identity server to grab the cookie
                // and write it to the .gstoken file.  This is required for the GS single-sign on model
                if (GSUtils.isGenomeSpace(url)) {
                    try {
                        java.util.List<HttpCookie> cookies = ((CookieManager) CookieManager.getDefault()).getCookieStore().get(url.toURI());
                        if (cookies != null) {
                            for (HttpCookie cookie : cookies) {
                                if (cookie.getName().equals("gs-token")) {
                                    GSUtils.setGSToken(cookie.getValue());
                                } else if (cookie.getName().equals("gs-username")) {
                                    GSUtils.setGSUser(cookie.getValue());
                                }
                            }
                        }
                    } catch (URISyntaxException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }

            // Redirects.  These can occur even if followRedirects == true if there is a change in protocol,
            // for example http -> https.
            if (code >= 300 && code < 400) {

                if (redirectCount > MAX_REDIRECTS) {
                    throw new IOException("Too many redirects");
                }

                String newLocation = conn.getHeaderField("Location");
                log.debug("Redirecting to " + newLocation);

                return openConnection(new URL(newLocation), requestProperties, method, redirectCount++);
            }

            // TODO -- handle other response codes.
            else if (code >= 400) {

                String message;
                if (code == 404) {
                    message = "File not found: " + url.toString();
                    throw new FileNotFoundException(message);
                } else {
                    message = conn.getResponseMessage();
                }
                String details = readErrorStream(conn);
                log.debug("error stream: " + details);
                log.debug(message);
                HttpResponseException exc = new HttpResponseException(code);

                throw exc;
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

    /**
     * The default authenticator
     */
    public class IGVAuthenticator extends Authenticator {

        /**
         * Called when password authentication is needed.
         *
         * @return
         */
        @Override
        protected PasswordAuthentication getPasswordAuthentication() {

            RequestorType type = getRequestorType();
            URL url = this.getRequestingURL();

            boolean isProxyChallenge = type == RequestorType.PROXY;
            if (isProxyChallenge) {
                if (proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                    return new PasswordAuthentication(proxySettings.user, proxySettings.pw.toCharArray());
                }
            }

            Frame owner = IGV.hasInstance() ? IGV.getMainFrame() : null;

            boolean isGenomeSpace = GSUtils.isGenomeSpace(url);

            LoginDialog dlg = new LoginDialog(owner, isGenomeSpace, url.toString(), isProxyChallenge);
            dlg.setVisible(true);
            if (dlg.isCanceled()) {
                return null;
            } else {
                final String userString = dlg.getUsername();
                final char[] userPass = dlg.getPassword();

                if (isProxyChallenge) {
                    proxySettings.user = userString;
                    proxySettings.pw = new String(userPass);
                }

                return new PasswordAuthentication(userString, userPass);
            }
        }
    }


    /**
     * Provide override for unit tests
     */
    public void setAuthenticator(Authenticator authenticator) {
        Authenticator.setDefault(authenticator);
    }

    /**
     * For unit tests
     */
    public void resetAuthenticator() {
        Authenticator.setDefault(new IGVAuthenticator());

    }


}
