/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util;

import biz.source_code.base64Coder.Base64Coder;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.ftp.FTPClient;
import htsjdk.samtools.util.ftp.FTPStream;
import org.apache.log4j.Logger;
import org.apache.tomcat.util.HttpDate;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.google.GoogleUtils;
import org.broad.igv.google.OAuthUtils;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.collections.CI;
import org.broad.igv.util.ftp.FTPUtils;

import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.*;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.time.ZonedDateTime;
import java.time.format.DateTimeParseException;
import java.util.List;
import java.util.*;
import java.util.zip.GZIPInputStream;

import static org.broad.igv.prefs.Constants.*;
import static org.broad.igv.util.stream.SeekableServiceStream.WEBSERVICE_URL;

/**
 * Wrapper utility class... for interacting with HttpURLConnection.
 *
 * @author Jim Robinson
 * @date 9/22/11
 */
public class HttpUtils {

    private static Logger log = Logger.getLogger(HttpUtils.class);

    private static HttpUtils instance;

    private Map<String, Boolean> byteRangeTestMap;

    private ProxySettings proxySettings = null;
    private final int MAX_REDIRECTS = 5;

    private String defaultUserName = null;
    private char[] defaultPassword = null;

    private Map<String, Collection<String>> headerMap = new HashMap<>();

    // static provided to support unit testing
    private static boolean BYTE_RANGE_DISABLED = false;
    private Map<URL, Boolean> headURLCache = new HashMap<URL, Boolean>();

    private class CachedRedirect {
        private URL url = null;
        private ZonedDateTime expires = null;
    }

    // remember HTTP redirects
    private final int DEFAULT_REDIRECT_EXPIRATION_MIN = 15;
    private Map<URL, CachedRedirect> redirectCache = new HashMap<URL, CachedRedirect>();

    /**
     * @return the single instance
     */
    public static HttpUtils getInstance() {
        if (instance == null) {
            instance = new HttpUtils();
        }
        return instance;
    }

    private HttpUtils() {

        disableCertificateValidation();
        Authenticator.setDefault(new IGVAuthenticator());

        try {
            System.setProperty("java.net.useSystemProxies", "true");
        } catch (Exception e) {
            log.info("Couldn't set useSystemProxies=true");
        }

        byteRangeTestMap = Collections.synchronizedMap(new HashMap());
    }

    /**
     * Create a URL from the given string.  Performs various mappings for google buckets,  amazon cNames, and
     * http -> https conversions
     *
     * @param urlString
     * @return
     * @throws MalformedURLException
     */
    public static URL createURL(String urlString) throws MalformedURLException {

        urlString = mapURL(urlString.trim());

        return new URL(urlString);
    }

    public static String mapURL(String urlString) throws MalformedURLException {

        if (urlString.startsWith("gs://")) {
            urlString = GoogleUtils.translateGoogleCloudURL(urlString);
        } else if (AmazonUtils.isAwsS3Path(urlString)) {
            try {
                urlString = AmazonUtils.translateAmazonCloudURL(urlString);
            } catch (IOException e) {
                log.error(e);
            }
        }

        if (GoogleUtils.isGoogleCloud(urlString)) {
            if (urlString.indexOf("alt=media") < 0) {
                urlString = URLUtils.addParameter(urlString, "alt=media");
            }
        }

        String host = URLUtils.getHost(urlString);
        if (host.equals("igv.broadinstitute.org")) {
            urlString = urlString.replace("igv.broadinstitute.org", "s3.amazonaws.com/igv.broadinstitute.org");
        } else if (host.equals("igvdata.broadinstitute.org")) {
            urlString = urlString.replace("igvdata.broadinstitute.org", "s3.amazonaws.com/igv.broadinstitute.org");
        } else if (host.equals("dn7ywbm9isq8j.cloudfront.net")) {
            urlString = urlString.replace("dn7ywbm9isq8j.cloudfront.net", "s3.amazonaws.com/igv.broadinstitute.org");
        } else if (host.equals("www.broadinstitute.org")) {
            urlString = urlString.replace("www.broadinstitute.org/igvdata", "data.broadinstitute.org/igvdata");
        } else if (host.equals("www.dropbox.com")) {
            urlString = urlString.replace("//www.dropbox.com", "//dl.dropboxusercontent.com");
        } else if (host.equals("drive.google.com")) {
            urlString = GoogleUtils.driveDownloadURL(urlString);
        }

        // data.broadinstitute.org requires https
        urlString = urlString.replace("http://data.broadinstitute.org", "https://data.broadinstitute.org");

        return urlString;
    }

    public static boolean isRemoteURL(String string) {
        return FileUtils.isRemote(string);
    }

    /**
     * Provided to support unit testing (force disable byte range requests)
     *
     * @return
     */
    public static void disableByteRange(boolean b) {
        BYTE_RANGE_DISABLED = b;
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
        return getContentsAsString(url, null);
    }


    public String getContentsAsString(URL url, Map<String, String> headers) throws IOException {

        InputStream is = null;

        HttpURLConnection conn = openConnection(url, headers);
        try {
            is = conn.getInputStream();
            return readContents(is);
        } catch (IOException e) {
            readErrorStream(conn);  // Consume content
            throw e;
        } finally {
            if (is != null) is.close();
        }
    }

    public String getContentsAsGzippedString(URL url) throws IOException {
        InputStream is = null;
        HttpURLConnection conn = openConnection(url, null);
        try {
            is = conn.getInputStream();
            return readContents(new GZIPInputStream(is));
        } catch (IOException e) {
            readErrorStream(conn);  // Consume content
            throw e;
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

        } catch (IOException e) {
            readErrorStream(conn);  // Consume content
            throw e;
        } finally {
            if (is != null) is.close();
        }
    }

    public String doPost(URL url, Map<String, String> params) throws IOException {

        StringBuilder postData = new StringBuilder();

        for (Map.Entry<String, String> param : params.entrySet()) {
            if (postData.length() != 0) postData.append('&');
            postData.append(param.getKey());
            postData.append('=');
            postData.append(param.getValue());
        }
        byte[] postDataBytes = postData.toString().getBytes();

        log.debug("Raw POST request: " + postData.toString());

        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("POST");
        conn.setRequestProperty("Content-Type", "application/x-www-form-urlencoded");
        conn.setDoOutput(true);
        conn.getOutputStream().write(postDataBytes);

        StringBuilder response = new StringBuilder();
        Reader in = new BufferedReader(new InputStreamReader(conn.getInputStream(), "UTF-8"));
        for (int c; (c = in.read()) >= 0; ) {
            response.append((char) c);
        }

        in.close();

        return response.toString();

    }

    /**
     * Open a connection stream for the URL.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public InputStream openConnectionStream(URL url) throws IOException {
        log.debug("Opening connection stream to  " + url);
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

        if (conn == null) {
            return null;
        }

        if ((requestProperties != null) && requestProperties.containsKey("Range") && conn.getResponseCode() != 206) {
            String msg = "Warning: range requested, but response code = " + conn.getResponseCode();
            log.error(msg);
        }

        try {
            InputStream input = conn.getInputStream();
            return input;
        } catch (IOException e) {
            readErrorStream(conn);  // Consume content
            throw e;
        }
    }

    public boolean resourceAvailable(String urlString) {

        URL url = null;
        try {
            url = HttpUtils.createURL(urlString);
        } catch (MalformedURLException e) {
            return false;
        }

        log.debug("Checking if resource is available: " + url);
        if (url.getProtocol().toLowerCase().equals("ftp")) {
            return FTPUtils.resourceAvailable(url);
        } else {
            HttpURLConnection conn = null;
            try {
                conn = openConnectionHeadOrGet(url);
                int code = conn.getResponseCode();
                return code >= 200 && code < 300;
            } catch (Exception e) {
                if (conn != null)
                    try {
                        readErrorStream(conn);  // Consume content
                    } catch (IOException e1) {
                        e1.printStackTrace();
                    }
                return false;
            } finally {
                if (conn != null) {
                    try {
                        conn.disconnect();
                    } catch (Exception e) {
                    }
                }
            }
        }
    }

    /**
     * First tries a HEAD request, then a GET request if the HEAD fails.
     * If the GET fails, the exception is thrown
     *
     * @param url
     * @return
     * @throws IOException
     */
    private HttpURLConnection openConnectionHeadOrGet(URL url) throws IOException {

        // Keep track of urls for which "HEAD" does not work (e.g. Amazon with signed urls).
        String urlString = url.toString();
        boolean isAWS = urlString.contains("AWSAccessKeyId");
        boolean tryHead =
                isAWS == false && (headURLCache.containsKey(url) ? headURLCache.get(url) : true);

        if (tryHead) {
            try {
                HttpURLConnection conn = openConnection(url, null, "HEAD");
                headURLCache.put(url, true);
                return conn;
            } catch (IOException e) {
                if (e instanceof FileNotFoundException) {
                    throw e;
                }
                log.debug("HEAD request failed for url: " + url.toExternalForm());
                log.debug("Trying GET instead for url: " + url.toExternalForm());
                headURLCache.put(url, false);
            }
        }
        return openConnection(url, null, "GET");
    }

    public String getHeaderField(URL url, String key) throws IOException {
        HttpURLConnection conn = openConnectionHeadOrGet(url);
        if (conn == null) return null;
        return conn.getHeaderField(key);
    }

    public long getLastModified(URL url) throws IOException {
        HttpURLConnection conn = openConnectionHeadOrGet(url);
        if (conn == null) return 0;
        return conn.getLastModified();
    }

    public long getContentLength(URL url) throws IOException {

        try {
            String contentLengthString = getHeaderField(url, "Content-Length");
            if (contentLengthString == null) {
                return -1;
            } else {
                return Long.parseLong(contentLengthString);
            }
        } catch (Exception e) {
            log.error("Error fetching content length", e);
            return -1;
        }
    }

    /**
     * Compare a local and remote resource, returning true if it is believed that the
     * remote file is newer than the local file
     *
     * @param file
     * @param url
     * @return true if the files are the same or the local file is newer, false if the remote file has been modified wrt the local one.
     * @throws IOException
     */
    public boolean remoteIsNewer(File file, URL url) throws IOException {

        if (!file.exists()) {
            return false;
        }

        HttpURLConnection conn = openConnection(url, null, "HEAD");

        // Compare last-modified dates
        String lastModifiedString = conn.getHeaderField("Last-Modified");
        if (lastModifiedString == null) {
            return true;                    // Assume its changed
        } else {
            HttpDate date = new HttpDate();
            date.parse(lastModifiedString);
            long remoteModifiedTime = date.getTime();
            long localModifiedTime = file.lastModified();
            return remoteModifiedTime > localModifiedTime;
        }


    }


    public void updateProxySettings() {
        boolean useProxy;
        String proxyHost;
        int proxyPort = -1;
        boolean auth = false;
        String user = null;
        String pw = null;

        IGVPreferences prefMgr = PreferencesManager.getPreferences();
        proxyHost = prefMgr.get(PROXY_HOST, null);
        try {
            proxyPort = Integer.parseInt(prefMgr.get(PROXY_PORT, "-1"));
        } catch (NumberFormatException e) {
            proxyPort = -1;
        }
        useProxy = prefMgr.getAsBoolean(USE_PROXY) && proxyHost != null && proxyHost.trim().length() > 0;
        auth = prefMgr.getAsBoolean(PROXY_AUTHENTICATE);
        user = prefMgr.get(PROXY_USER, null);
        String pwString = prefMgr.get(PROXY_PW, null);
        if (pwString != null) {
            pw = Utilities.base64Decode(pwString);
        }

        String proxyTypeString = prefMgr.get(PROXY_TYPE, "HTTP");
        Proxy.Type type = Proxy.Type.valueOf(proxyTypeString.trim().toUpperCase());

        String proxyWhitelistString = prefMgr.get(PROXY_WHITELIST);
        Set<String> whitelist = proxyWhitelistString == null ? new HashSet<String>() :
                new HashSet(Arrays.asList(Globals.commaPattern.split(proxyWhitelistString)));

        proxySettings = new ProxySettings(useProxy, user, pw, auth, proxyHost, proxyPort, type, whitelist);
    }

    /**
     * Get the system defined proxy defined for the URI, or null if
     * not available. May also return a {@code Proxy} object which
     * represents a direct connection
     *
     * @param uri
     * @return
     */
    private Proxy getSystemProxy(String uri) {
        try {
            if (PreferencesManager.getPreferences().getAsBoolean("DEBUG.PROXY"))
                log.info("Getting system proxy for " + uri);
            ProxySelector selector = ProxySelector.getDefault();
            List<Proxy> proxyList = selector.select(new URI(uri));
            return proxyList.get(0);
        } catch (URISyntaxException e) {
            log.error(e.getMessage(), e);
            return null;
        } catch (NullPointerException e) {
            return null;
        } catch (Exception e) {
            log.error(e.getMessage(), e);
            return null;
        }

    }


    /**
     * Calls {@link #downloadFile(String, java.io.File, Frame, String)}
     * with {@code dialogsParent = null, title = null}
     *
     * @param url
     * @param outputFile
     * @return RunnableResult
     * @throws IOException
     */
    public RunnableResult downloadFile(String url, File outputFile) throws IOException {
        URLDownloader downloader = downloadFile(url, outputFile, null, null);
        return downloader.getResult();
    }

    /**
     * @param url
     * @param outputFile
     * @param dialogsParent Parent of dialog to show progress. If null, none shown
     * @return URLDownloader used to perform download
     * @throws IOException
     */
    public URLDownloader downloadFile(String url, File outputFile, Frame dialogsParent, String dialogTitle) throws IOException {
        final URLDownloader urlDownloader = new URLDownloader(url, outputFile);
        boolean showProgressDialog = dialogsParent != null;
        if (!showProgressDialog) {
            urlDownloader.run();
            return urlDownloader;
        } else {
            javax.swing.ProgressMonitor monitor = new javax.swing.ProgressMonitor(IGV.getInstance().getMainPanel(), "Downloading " + outputFile.getName(), "", 0, 100);
            urlDownloader.setMonitor(monitor);
            ActionListener buttonListener = new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    urlDownloader.cancel(true);
                }
            };
            //  String permText = "Downloading " + url;
            //  String title = dialogTitle != null ? dialogTitle : permText;
            //  CancellableProgressDialog dialog = CancellableProgressDialog.showCancellableProgressDialog(dialogsParent, title, buttonListener, false, monitor);
            //  dialog.setPermText(permText);

            //  Dimension dms = new Dimension(600, 150);
            //  dialog.setPreferredSize(dms);
            //  dialog.setSize(dms);
            //  dialog.validate();

            LongRunningTask.submit(urlDownloader);
            return urlDownloader;
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
        return ParsingUtils.readContentsFromStream(is);
    }

    public String readErrorStream(HttpURLConnection connection) throws IOException {
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

    public HttpURLConnection delete(URL url) throws IOException {
        return openConnection(url, Collections.<String, String>emptyMap(), "DELETE");
    }

    public HttpURLConnection openConnection(URL url, Map<String, String> requestProperties) throws IOException {
        return openConnection(url, requestProperties, "GET");
    }

    private HttpURLConnection openConnection(URL url, Map<String, String> requestProperties, String method) throws IOException {
        return openConnection(url, requestProperties, method, 0, 0);
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

            URL url, Map<String, String> requestProperties, String method, int redirectCount, int retries) throws IOException {

        // if we're already seen a redirect for this URL, use the updated one
        if (redirectCache.containsKey(url)) {
            CachedRedirect cr = redirectCache.get(url);
            if (ZonedDateTime.now().compareTo(cr.expires) < 0.0) {
                // now() is before our expiration
                log.debug("Found URL in redirection cache: " + url + " ->" + redirectCache.get(url).url);
                url = cr.url;
            } else {
                log.debug("Removing expired URL from redirection cache: " + url);
                redirectCache.remove(url);
            }

        }

        // if the url points to a openid location instead of a oauth2.0 location, used the fina and replace
        // string to dynamically map url - dwm08
        if (url.getHost().equals(GoogleUtils.GOOGLE_API_HOST) && OAuthUtils.findString != null && OAuthUtils.replaceString != null) {
            url = HttpUtils.createURL(url.toExternalForm().replaceFirst(OAuthUtils.findString, OAuthUtils.replaceString));
        }

        //Encode query string portions
        url = StringUtils.encodeURLQueryString(url);
        if (log.isTraceEnabled()) {
            log.trace(url);
        }

        //Encode base portions. Right now just spaces, most common case
        //TODO This is a hack and doesn't work for all characters which need it
        if (StringUtils.countChar(url.toExternalForm(), ' ') > 0) {
            String newPath = url.toExternalForm().replaceAll(" ", "%20");
            url = HttpUtils.createURL(newPath);
        }

        Proxy sysProxy = null;
        boolean igvProxySettingsExist = proxySettings != null && proxySettings.useProxy;
        boolean checkSystemProxy =
                !PreferencesManager.getPreferences().getAsBoolean("PROXY.DISABLE_CHECK") && !igvProxySettingsExist;

        //Only check for system proxy if igv proxy settings not found
        if (checkSystemProxy) {
            sysProxy = getSystemProxy(url.toExternalForm());
        }

        boolean useProxy =
                (sysProxy != null && sysProxy.type() != Proxy.Type.DIRECT) ||
                        (igvProxySettingsExist && !proxySettings.getWhitelist().contains(url.getHost()));

        HttpURLConnection conn;
        if (useProxy) {
            Proxy proxy = sysProxy;
            if (igvProxySettingsExist) {
                if (proxySettings.type == Proxy.Type.DIRECT) {

                    if (PreferencesManager.getPreferences().getAsBoolean("DEBUG.PROXY")) {
                        log.info("NO_PROXY");
                    }

                    proxy = Proxy.NO_PROXY;
                } else {
                    if (PreferencesManager.getPreferences().getAsBoolean("DEBUG.PROXY")) {
                        log.info("PROXY " + proxySettings.proxyHost + "  " + proxySettings.proxyPort);
                    }

                    proxy = new Proxy(proxySettings.type, new InetSocketAddress(proxySettings.proxyHost, proxySettings.proxyPort));
                }
            }
            conn = (HttpURLConnection) url.openConnection(proxy);
            if (igvProxySettingsExist && proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                byte[] bytes = (proxySettings.user + ":" + proxySettings.pw).getBytes();

                String encodedUserPwd = String.valueOf(Base64Coder.encode(bytes));
                conn.setRequestProperty("Proxy-Authorization", "Basic " + encodedUserPwd);
            }
        } else {
            if (PreferencesManager.getPreferences().getAsBoolean("DEBUG.PROXY")) {
                log.info("PROXY NOT USED ");
                if (proxySettings.getWhitelist().contains(url.getHost())) {
                    log.info(url.getHost() + " is whitelisted");
                }
                ;
            }
            conn = (HttpURLConnection) url.openConnection();
        }


        if (!"HEAD".equals(method)) {
            conn.setRequestProperty("Accept", "text/plain");
        }

        conn.setConnectTimeout(Globals.CONNECT_TIMEOUT);
        conn.setReadTimeout(Globals.READ_TIMEOUT);
        conn.setRequestMethod(method);
        conn.setRequestProperty("Connection", "Keep-Alive");
        // we'll handle redirects manually, allowing us to cache the new URL
        conn.setInstanceFollowRedirects(false);

        if (requestProperties != null) {
            for (Map.Entry<String, String> prop : requestProperties.entrySet()) {
                conn.setRequestProperty(prop.getKey(), prop.getValue());
            }
        }

        Collection<String> headers = headerMap.get(url.getHost());
        if (headers != null) {
            for (String h : headers) {
                String[] kv = h.split(":");
                if (kv.length == 2) {
                    conn.setRequestProperty(kv[0], kv[1]);
                }
            }
        }

        conn.setRequestProperty("User-Agent", Globals.applicationString());

        // If this is a Google URL and we have an access token use it.
        if (GoogleUtils.isGoogleURL(url.toExternalForm())) {
            String token = OAuthUtils.getInstance().getProvider().getAccessToken();
            if (token != null) {
                conn.setRequestProperty("Authorization", "Bearer " + token);
            }
            if (GoogleUtils.getProjectID() != null &&
                    GoogleUtils.getProjectID().length() > 0 &&
                    !hasQueryParameter(url, "userProject")) {
                url = addQueryParameter(url, "userProject", GoogleUtils.getProjectID());
            }
        }

        if (method.equals("PUT")) {
            return conn;
        } else {

            int code = conn.getResponseCode();

            if (!isDropboxHost(url.getHost()) && requestProperties != null && requestProperties.containsKey("Range") && code == 200 && method.equals("GET")) {

                log.error("Range header removed by client or ignored by server for url: " + url.toString());

                if (!SwingUtilities.isEventDispatchThread()) {
                    MessageUtils.showMessage("Warning: unsuccessful attempt to execute 'Range byte' request to host " + url.getHost());
                }

                byteRangeTestMap.put(url.getHost(), false);
                String[] positionString = requestProperties.get("Range").split("=")[1].split("-");
                int length = Integer.parseInt(positionString[1]) - Integer.parseInt(positionString[0]) + 1;
                requestProperties.remove("Range"); // < VERY IMPORTANT
                URL wsUrl = HttpUtils.createURL(WEBSERVICE_URL + "?file=" + url.toExternalForm() + "&position=" + positionString[0] + "&length=" + length);
                return openConnection(wsUrl, requestProperties, "GET", redirectCount, retries);
            }

            if (log.isDebugEnabled()) {
                //logHeaders(conn);
            }

            // Redirects.  These can occur even if followRedirects == true if there is a change in protocol,
            // for example http -> https.
            if (code >= 300 && code < 400) {
                if (redirectCount > MAX_REDIRECTS) {
                    throw new IOException("Too many redirects");
                }

                CachedRedirect cr = new CachedRedirect();
                cr.url = new URL(conn.getHeaderField("Location"));
                if (cr.url != null) {
                    cr.expires = ZonedDateTime.now().plusMinutes(DEFAULT_REDIRECT_EXPIRATION_MIN);
                    String s;
                    if ((s = conn.getHeaderField("Cache-Control")) != null) {

                        // cache-control takes priority
                        CacheControl cc = null;
                        try {
                            cc = CacheControl.valueOf(s);
                        } catch (IllegalArgumentException e) {
                            // use default
                        }
                        if (cc != null) {
                            if (cc.isNoCache()) {
                                // set expires to null, preventing caching
                                cr.expires = null;
                            } else if (cc.getMaxAge() > 0) {
                                cr.expires = ZonedDateTime.now().plusSeconds(cc.getMaxAge());
                            }
                        }
                    } else if ((s = conn.getHeaderField("Expires")) != null) {
                        // no cache-control header, so try "expires" next
                        try {
                            cr.expires = ZonedDateTime.parse(s);
                        } catch (DateTimeParseException e) {
                            // use default
                        }
                    }
                    if (cr.expires != null) {
                        redirectCache.put(url, cr);
                        log.debug("Redirecting to " + cr.url);
                        return openConnection(HttpUtils.createURL(cr.url.toString()), requestProperties, method, ++redirectCount, retries);
                    }
                }
            }

            // TODO -- handle other response codes.
            else if (code >= 400) {

                String message;

                // TODO -- detect Google requestor pay failure

                if (code == 404) {
                    message = "File not found: " + url.toString();
                    throw new FileNotFoundException(message);
                } else if (code == 401) {
                    if (GoogleUtils.isGoogleURL(url.toExternalForm()) && retries == 0) {
                        GoogleUtils.checkLogin();
                        return openConnection(url, requestProperties, method, redirectCount, ++retries);
                    }
                    message = "You must log in to access this file";
                    throw new HttpResponseException(code, message, "");
                } else if (code == 403) {
                    message = "Access forbidden";
                    throw new HttpResponseException(code, message, "");
                } else if (code == 416) {
                    throw new UnsatisfiableRangeException(conn.getResponseMessage());
                } else {
                    message = conn.getResponseMessage();
                    String details = readErrorStream(conn);

                    if (url.getHost().equals("www.googleapis.com") && details.contains("requester pays bucket")) {
                        MessageUtils.showMessage("<html>" + details + "<br>Use Google menu to set project.");
                    }

                    throw new HttpResponseException(code, message, details);
                }
            }
        }
        return conn;
    }

    private boolean isDropboxHost(String host) {
        return (host.equals("dl.dropboxusercontent.com") || host.equals("www.dropbox.com"));
    }

    private URL addQueryParameter(URL url, String userProject, String projectID) {
        String urlString = url.toExternalForm();
        urlString = urlString + (urlString.contains("?") ? "&" : "?") + userProject + "=" + projectID;
        try {
            return new URL(urlString);
        } catch (MalformedURLException e) {
            log.error("Error adding query parameter", e);
            return url;
        }
    }

    private boolean hasQueryParameter(URL url, String parameter) {
        String urlstring = url.toExternalForm();
        if (urlstring.contains("?")) {
            int idx = urlstring.indexOf('?');
            return urlstring.substring(idx).contains("parameter" + "=");
        } else {
            return false;
        }
    }


    //Used for testing sometimes, please do not delete
    private void logHeaders(HttpURLConnection conn) {
        Map<String, List<String>> headerFields = conn.getHeaderFields();
        log.debug("Headers for " + conn.getURL());
        for (Map.Entry<String, List<String>> header : headerFields.entrySet()) {
            log.debug(header.getKey() + ": " + StringUtils.join(header.getValue(), ","));
        }
    }


    public void setDefaultPassword(String defaultPassword) {
        this.defaultPassword = defaultPassword.toCharArray();
    }

    public void setDefaultUserName(String defaultUserName) {
        this.defaultUserName = defaultUserName;
    }

    public void clearDefaultCredentials() {
        this.defaultPassword = null;
        this.defaultUserName = null;
    }


    /**
     * Test to see if this client can successfully retrieve a portion of a remote file using the byte-range header.
     * This is not a test of the server, but the client.  In some environments the byte-range header gets removed
     * by filters after the request is made by IGV.
     *
     * @return
     */
    public boolean useByteRange(URL url) throws IOException {

        if (BYTE_RANGE_DISABLED) return false;

        // We can test byte-range success for hosts we can reach.
        synchronized (byteRangeTestMap) {
            final String host = url.getHost();

            if (byteRangeTestMap.containsKey(host)) {
                return byteRangeTestMap.get(host);
            } else {
                SeekableStream str = null;
                try {

                    boolean byteRangeTestSuccess = testByteRange(url);

                    if (byteRangeTestSuccess) {
                        log.info("Range-byte request succeeded");
                    } else {
                        log.info("Range-byte test failed -- Host: " + host +
                                " does not support range-byte requests or there is a problem with client network environment.");
                    }

                    byteRangeTestMap.put(host, byteRangeTestSuccess);
                    return byteRangeTestSuccess;


                } finally {
                    if (str != null) try {
                        str.close();
                    } catch (IOException e) {
                        log.error("Error closing stream (" + url.toExternalForm() + ")", e);
                    }
                }
            }

        }
    }

    public boolean testByteRange(URL url) throws IOException {
        Map<String, String> params = new HashMap();
        String byteRange = "bytes=" + 0 + "-" + 10;
        params.put("Range", byteRange);
        HttpURLConnection conn = HttpUtils.getInstance().openConnection(url, params);
        int statusCode = conn.getResponseCode();
        boolean byteRangeTestSuccess = (statusCode == 206);
        readFully(conn.getInputStream(), new byte[10]);
        return byteRangeTestSuccess;
    }

    /**
     * Add an http header string to be applied the the specified URLs.  Used to support command line specification
     * of authentication headers
     *
     * @param headers
     * @param urls
     */
    public void addHeaders(Collection<String> headers, List<String> urls) {
        for (String u : urls) {
            if (isRemoteURL(u)) {
                try {
                    URL url = new URL(mapURL(u));
                    headerMap.put(url.getHost(), headers);
                    System.out.println("Added " + url.getHost() + " -> " + headers);

                } catch (MalformedURLException e) {
                    log.error("Error parsing URL " + u, e);
                }
            }
        }
    }


    private String stripParameters(String url) {
        int idx = url.indexOf("?");
        if (idx > 0) {
            return url.substring(0, idx);
        } else {
            return url;
        }
    }

    public void shutdown() {
        // Do any cleanup required here
    }

    public static class ProxySettings {
        boolean auth = false;
        String user;
        String pw;
        boolean useProxy;
        String proxyHost;
        int proxyPort = -1;
        Proxy.Type type;
        Set<String> whitelist;


        public ProxySettings(boolean useProxy, String user, String pw, boolean auth, String proxyHost, int proxyPort,
                             Proxy.Type proxyType, Set<String> whitelist) {
            this.auth = auth;
            this.proxyHost = proxyHost;
            this.proxyPort = proxyPort;
            this.pw = pw;
            this.useProxy = useProxy;
            this.user = user;
            this.type = proxyType;
            this.whitelist = whitelist;
        }

        public Set<String> getWhitelist() {
            return whitelist;
        }
    }

    /**
     * The default authenticator
     */
    public class IGVAuthenticator extends Authenticator {

        Hashtable<String, PasswordAuthentication> pwCache = new Hashtable<String, PasswordAuthentication>();
        HashSet<String> cacheAttempts = new HashSet<String>();

        /**
         * Called when password authentication is needed.
         *
         * @return
         */
        @Override
        protected synchronized PasswordAuthentication getPasswordAuthentication() {

            RequestorType type = getRequestorType();
            String urlString = getRequestingURL().toString();
            boolean isProxyChallenge = type == RequestorType.PROXY;

            // Cache user entered PWs.  In normal use this shouldn't be necessary as credentials are cached upstream,
            // but if loading many files in parallel (e.g. from sessions) calls to this method can queue up before the
            // user enters their credentials, causing needless reentry.
            String pKey = type.toString() + getRequestingProtocol() + getRequestingHost();
            PasswordAuthentication pw = pwCache.get(pKey);
            if (pw != null) {
                // Prevents infinite loop if credentials are incorrect
                if (cacheAttempts.contains(urlString)) {
                    cacheAttempts.remove(urlString);
                } else {
                    cacheAttempts.add(urlString);
                    return pw;
                }
            }

            if (isProxyChallenge) {
                if (proxySettings.auth && proxySettings.user != null && proxySettings.pw != null) {
                    return new PasswordAuthentication(proxySettings.user, proxySettings.pw.toCharArray());
                }
            }

            if (defaultUserName != null && defaultPassword != null) {
                return new PasswordAuthentication(defaultUserName, defaultPassword);
            }

            Frame owner = IGV.hasInstance() ? IGV.getMainFrame() : null;

            LoginDialog dlg = new LoginDialog(owner, urlString, isProxyChallenge);
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

                pw = new PasswordAuthentication(userString, userPass);
                pwCache.put(pKey, pw);
                return pw;
            }
        }
    }


    static boolean isExpectedRangeMissing(URLConnection conn, Map<String, String> requestProperties) {
        final boolean rangeRequested = (requestProperties != null) && (new CI.CIHashMap<String>(requestProperties)).containsKey("Range");
        if (!rangeRequested) return false;

        Map<String, List<String>> headerFields = conn.getHeaderFields();
        boolean rangeReceived = (headerFields != null) && (new CI.CIHashMap<List<String>>(headerFields)).containsKey("Content-Range");
        return !rangeReceived;
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


    /**
     * Useful helper function
     */
    public static void readFully(InputStream is, byte b[]) throws IOException {
        int len = b.length;
        if (len < 0) {
            throw new IndexOutOfBoundsException();
        }
        int n = 0;
        while (n < len) {
            int count = is.read(b, n, len - n);
            if (count < 0) {
                throw new EOFException();
            }
            n += count;
        }
    }

    public class UnsatisfiableRangeException extends RuntimeException {

        String message;

        public UnsatisfiableRangeException(String message) {
            super(message);
            this.message = message;
        }
    }

    static class CacheControl {

        boolean noCache = false;
        long maxAge = 0;

        static CacheControl valueOf(String s) {
            CacheControl cc = new CacheControl();
            String[] tokens = Globals.commaPattern.split(s);
            for (String t : tokens) {
                t = t.trim().toLowerCase();
                if (t.startsWith("no-cache")) {
                    cc.noCache = true;
                } else if (t.startsWith("max-age")) {
                    String[] ma = Globals.equalPattern.split(t);
                    cc.maxAge = Long.parseLong(ma[1].trim());
                }
            }
            return cc;
        }

        private CacheControl() {
        }

        public boolean isNoCache() {
            return noCache;
        }

        public long getMaxAge() {
            return maxAge;
        }
    }
}
