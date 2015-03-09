/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.util;

import biz.source_code.base64Coder.Base64Coder;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.ftp.FTPClient;
import htsjdk.samtools.util.ftp.FTPStream;
import org.apache.log4j.Logger;
import org.apache.tomcat.util.HttpDate;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.CancellableProgressDialog;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.collections.CI;
import org.broad.igv.util.ftp.FTPUtils;
import org.broad.igv.util.stream.IGVSeekableHTTPStream;
import org.broad.igv.util.stream.IGVUrlHelper;

import javax.net.ssl.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.*;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

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
    private static Pattern URLmatcher = Pattern.compile(".{1,8}://.*");

    // static provided to support unit testing
    private static boolean BYTE_RANGE_DISABLED = false;
    private Map<URL, Boolean> headURLCache = new HashMap<URL, Boolean>();

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

        htsjdk.tribble.util.ParsingUtils.registerHelperClass(IGVUrlHelper.class);

        disableCertificateValidation();
        CookieHandler.setDefault(new IGVCookieManager());
        Authenticator.setDefault(new IGVAuthenticator());

        try {
            System.setProperty("java.net.useSystemProxies", "true");
        } catch (Exception e) {
            log.info("Couldn't set useSystemProxies=true");
        }

        byteRangeTestMap = Collections.synchronizedMap(new HashMap());
    }

    public static boolean isRemoteURL(String string) {
        String lcString = string.toLowerCase();
        return lcString.startsWith("http://") || lcString.startsWith("https://") || lcString.startsWith("ftp://");
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

        InputStream is = null;
        HttpURLConnection conn = openConnection(url, null);
        try {
            is = getInputStream(conn);
            return readContents(is);

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
            is = getInputStream(conn);
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

        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("POST");
        conn.setRequestProperty("Content-Type", "application/x-www-form-urlencoded");
        //conn.setRequestProperty("Content-Length", String.valueOf(postDataBytes.length));
        conn.setDoOutput(true);
        conn.getOutputStream().write(postDataBytes);

        StringBuilder response = new StringBuilder();
        Reader in = new BufferedReader(new InputStreamReader(getInputStream(conn), "UTF-8"));
        for (int c; (c = in.read()) >= 0; ) {
            response.append((char) c);
        }
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

        boolean rangeRequestedNotReceived = isExpectedRangeMissing(conn, requestProperties);
        if (rangeRequestedNotReceived) {
            String msg = "Byte range requested, but no Content-Range header in response";
            log.error(msg);
//            if(Globals.isTesting()){
//                throw new IOException(msg);
//            }
        }

        try {
            InputStream input = getInputStream(conn);
            return input;
        } catch (IOException e) {
            readErrorStream(conn);  // Consume content
            throw e;
        }
    }

    private InputStream getInputStream(HttpURLConnection conn) throws IOException {
        InputStream input = conn.getInputStream();
        //  if ("gzip".equals(conn.getContentEncoding())) {
        //      input = new GZIPInputStream(input);
        //  }
        return input;
    }

    boolean isExpectedRangeMissing(URLConnection conn, Map<String, String> requestProperties) {
        final boolean rangeRequested = (requestProperties != null) && (new CI.CIHashMap<String>(requestProperties)).containsKey("Range");
        if (!rangeRequested) return false;

        Map<String, List<String>> headerFields = conn.getHeaderFields();
        boolean rangeReceived = (headerFields != null) && (new CI.CIHashMap<List<String>>(headerFields)).containsKey("Content-Range");
        return !rangeReceived;
    }

    public boolean resourceAvailable(URL url) {
        log.debug("Checking if resource is available: " + url);
        if (url.getProtocol().toLowerCase().equals("ftp")) {
            return FTPUtils.resourceAvailable(url);
        } else {
            try {
                HttpURLConnection conn = openConnectionHeadOrGet(url);
                int code = conn.getResponseCode();
                return code >= 200 && code < 300;
            } catch (Exception e) {
                return false;
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
        boolean tryHead = headURLCache.containsKey(url) ? headURLCache.get(url) : true;

        if (tryHead) {
            try {
                HttpURLConnection conn = openConnection(url, null, "HEAD");
                headURLCache.put(url, true);
                return conn;
            } catch (IOException e) {
                if (e instanceof FileNotFoundException) {
                    throw e;
                }
                log.info("HEAD request failed for url: " + url.toExternalForm() + ".  Trying GET");
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
     * @param compareContentLength Whether to use the content length to compare files. If false, only
     *                             the modified date is used
     * @return true if the files are the same or the local file is newer, false if the remote file has been modified wrt the local one.
     * @throws IOException
     */
    public boolean remoteIsNewer(File file, URL url, boolean compareContentLength) throws IOException {

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
            return true;
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

        String proxyTypeString = prefMgr.get(PreferenceManager.PROXY_TYPE, "HTTP");
        Proxy.Type type = Proxy.Type.valueOf(proxyTypeString.trim().toUpperCase());

        proxySettings = new ProxySettings(useProxy, user, pw, auth, proxyHost, proxyPort, type);
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
            log.debug("Getting system proxy for " + uri);
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
            ProgressMonitor monitor = new ProgressMonitor();
            urlDownloader.setMonitor(monitor);
            ActionListener buttonListener = new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    urlDownloader.cancel(true);
                }
            };
            String permText = "Downloading " + url;
            String title = dialogTitle != null ? dialogTitle : permText;
            CancellableProgressDialog dialog = CancellableProgressDialog.showCancellableProgressDialog(dialogsParent, title, buttonListener, false, monitor);
            dialog.setPermText(permText);

            Dimension dms = new Dimension(600, 150);
            dialog.setPreferredSize(dms);
            dialog.setSize(dms);
            dialog.validate();

            LongRunningTask.submit(urlDownloader);
            return urlDownloader;
        }
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
            inputStream = getInputStream(urlconnection);
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
                        return new java.security.cert.X509Certificate[0];
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

        // Install a callback to allow the igv Amazon and local hosts
        HttpsURLConnection.setDefaultHostnameVerifier(new HostnameVerifier() {
            public boolean verify(String hostname, SSLSession session) {
                if (hostname.equals("igv.broadinstitute.org") || hostname.equals("igvdata.broadinstitute.org") || hostname.equals("localhost"))
                    return true;
                return false;
            }
        });

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

    public HttpURLConnection delete(URL url) throws IOException {
        return openConnection(url, Collections.<String, String>emptyMap(), "DELETE");
    }

    HttpURLConnection openConnection(URL url, Map<String, String> requestProperties) throws IOException {
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

        //Encode query string portions
        url = StringUtils.encodeURLQueryString(url);
        if (log.isTraceEnabled()) {
            log.trace(url);
        }

        //Encode base portions. Right now just spaces, most common case
        //TODO This is a hack and doesn't work for all characters which need it
        if (StringUtils.countChar(url.toExternalForm(), ' ') > 0) {
            String newPath = url.toExternalForm().replaceAll(" ", "%20");
            url = new URL(newPath);
        }

        Proxy sysProxy = null;
        boolean igvProxySettingsExist = proxySettings != null && proxySettings.useProxy;
        //Only check for system proxy if igv proxy settings not found
        if (!igvProxySettingsExist) {
            sysProxy = getSystemProxy(url.toExternalForm());
        }
        boolean useProxy = sysProxy != null || igvProxySettingsExist;

        HttpURLConnection conn;
        if (useProxy) {
            Proxy proxy = sysProxy;
            if (igvProxySettingsExist) {
                if (proxySettings.type == Proxy.Type.DIRECT) {
                    proxy = Proxy.NO_PROXY;
                } else {
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
            conn = (HttpURLConnection) url.openConnection();
        }

        if (GSUtils.isGenomeSpace(url)) {
            conn.setRequestProperty("Accept", "application/json,text/plain");
        } else {
            if (!"HEAD".equals(method))
                conn.setRequestProperty("Accept", "text/plain");
        }

        //------//
        //There seems to be a bug with JWS caches
        //So we avoid caching

        //This default is persistent, really should be available statically but isn't
        conn.setDefaultUseCaches(false);
        conn.setUseCaches(false);
        //------//


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

            if (log.isDebugEnabled()) {
                //logHeaders(conn);
            }

            // Redirects.  These can occur even if followRedirects == true if there is a change in protocol,
            // for example http -> https.
            if (code >= 300 && code < 400) {

                if (redirectCount > MAX_REDIRECTS) {
                    throw new IOException("Too many redirects");
                }

                String newLocation = conn.getHeaderField("Location");
                log.debug("Redirecting to " + newLocation);

                return openConnection(new URL(newLocation), requestProperties, method, ++redirectCount);
            }

            // TODO -- handle other response codes.
            else if (code >= 400) {

                String message;
                if (code == 404) {
                    message = "File not found: " + url.toString();
                    throw new FileNotFoundException(message);
                } else if (code == 401) {
                    // Looks like this only happens when user hits "Cancel".
                    // message = "Not authorized to view this file";
                    // JOptionPane.showMessageDialog(null, message, "HTTP error", JOptionPane.ERROR_MESSAGE);
                    redirectCount = MAX_REDIRECTS + 1;
                    return null;
                } else {
                    message = conn.getResponseMessage();
                }
                String details = readErrorStream(conn);
                log.error("URL: " + url.toExternalForm() + ". error stream: " + details);
                log.error("Code: " + code + ". " + message);
                HttpResponseException exc = new HttpResponseException(code);
                throw exc;
            }
        }
        return conn;
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
    public boolean useByteRange(URL url) {

        if (BYTE_RANGE_DISABLED) return false;

        // We can test byte-range success for hosts we can reach.
        synchronized (byteRangeTestMap) {
            final String host = url.getHost();

            if (byteRangeTestMap.containsKey(host)) {
                return byteRangeTestMap.get(host);
            } else {
                SeekableStream str = null;
                try {
                    boolean byteRangeTestSuccess = true;

                    if (host.contains("broadinstitute.org")) {
                        byteRangeTestSuccess = testBroadHost(host);
                    } else {
                        // Non-broad URL
                        int l = (int) Math.min(1000, HttpUtils.getInstance().getContentLength(url));
                        if (l > 100) {

                            byte[] firstBytes = new byte[l];
                            str = new IGVSeekableHTTPStream(url);
                            str.readFully(firstBytes);

                            int end = firstBytes.length;
                            int start = end - 100;
                            str.seek(start);
                            int len = end - start;
                            byte[] buffer = new byte[len];
                            int n = 0;
                            while (n < len) {
                                int count = str.read(buffer, n, len - n);
                                if (count < 0)
                                    throw new EOFException();
                                n += count;
                            }

                            for (int i = 0; i < len; i++) {
                                if (buffer[i] != firstBytes[i + start]) {
                                    byteRangeTestSuccess = false;
                                    break;
                                }
                            }
                        } else {
                            // Too small a sample to test, or unknown content length.  Return "true" but don't record
                            // this host as tested.
                            return true;
                        }
                    }

                    if (byteRangeTestSuccess) {
                        log.info("Range-byte request succeeded");
                    } else {
                        log.info("Range-byte test failed -- problem with client network environment.");
                    }

                    byteRangeTestMap.put(host, byteRangeTestSuccess);
                    return byteRangeTestSuccess;


                } catch (IOException e) {
                    log.error("Error while testing byte range " + e.getMessage());
                    // We could not reach the test server, so we can't know if this client can do byte-range tests or
                    // not.  Take the "optimistic" view.
                    return true;
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

    private boolean testBroadHost(String host) throws IOException {
        // Test broad urls for successful byte range requests.
        log.info("Testing range-byte request on host: " + host);

        String testURL;
        if (host.startsWith("igvdata.broadinstitute.org")) {
            // Amazon cloud front
            testURL = "http://igvdata.broadinstitute.org/genomes/seq/hg19/chr12.txt";
        } else if (host.startsWith("igv.broadinstitute.org")) {
            // Amazon S3
            testURL = "http://igv.broadinstitute.org/genomes/seq/hg19/chr12.txt";
        } else {
            // All others
            testURL = "http://www.broadinstitute.org/igvdata/annotations/seq/hg19/chr12.txt";
        }

        byte[] expectedBytes = {'T', 'C', 'G', 'C', 'T', 'T', 'G', 'A', 'A', 'C', 'C', 'C', 'G', 'G',
                'G', 'A', 'G', 'A', 'G', 'G'};
        IGVSeekableHTTPStream str = null;

        try {
            str = new IGVSeekableHTTPStream(new URL(testURL));

            str.seek(25350000);
            byte[] buffer = new byte[80000];
            str.read(buffer);
            String result = new String(buffer);
            for (int i = 0; i < expectedBytes.length; i++) {
                if (buffer[i] != expectedBytes[i]) {
                    return false;
                }
            }
            return true;
        } finally {
            if (str != null) str.close();
        }
    }


    public void shutdown() {
        // Do any cleanup required here
    }

    /**
     * Checks if the string is a URL (not necessarily remote, can be any protocol)
     *
     * @param f
     * @return
     */
    public static boolean isURL(String f) {
        return f.startsWith("http:") || f.startsWith("ftp:") || f.startsWith("https:") || URLmatcher.matcher(f).matches();
    }

    public static Map<String, String> parseQueryString(String query) {
        String[] params = query.split("&");
        Map<String, String> map = new HashMap<String, String>();
        for (String param : params) {
            String[] name_val = param.split("=", 2);
            if (name_val.length == 2) {
                map.put(name_val[0], name_val[1]);
            }
        }
        return map;
    }


    public static class ProxySettings {
        boolean auth = false;
        String user;
        String pw;
        boolean useProxy;
        String proxyHost;
        int proxyPort = -1;
        Proxy.Type type;

        public ProxySettings(boolean useProxy, String user, String pw, boolean auth, String proxyHost, int proxyPort) {
            this(useProxy, user, pw, auth, proxyHost, proxyPort, Proxy.Type.HTTP);
            this.auth = auth;
        }

        public ProxySettings(boolean useProxy, String user, String pw, boolean auth, String proxyHost, int proxyPort, Proxy.Type proxyType) {
            this.auth = auth;
            this.proxyHost = proxyHost;
            this.proxyPort = proxyPort;
            this.pw = pw;
            this.useProxy = useProxy;
            this.user = user;
            this.type = proxyType;
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

            boolean isGenomeSpace = GSUtils.isGenomeSpace(getRequestingURL());
            if (isGenomeSpace) {
                // If we are being challenged by GS the token must be bad/expired
                GSUtils.logout();
            }

            LoginDialog dlg = new LoginDialog(owner, isGenomeSpace, urlString, isProxyChallenge);
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

    /**
     * Runnable for downloading a file from a URL.
     * Downloading is buffered, and can be cancelled (between buffers)
     * via {@link #cancel(boolean)}
     */
    public class URLDownloader implements Runnable {

        private ProgressMonitor monitor = null;

        private final URL srcUrl;
        private final File outputFile;

        private volatile boolean started = false;
        private volatile boolean done = false;
        private volatile boolean cancelled = false;
        private volatile RunnableResult result;

        public URLDownloader(String url, File outputFile) throws MalformedURLException {
            this.srcUrl = new URL(url);
            this.outputFile = outputFile;
        }

        @Override
        public void run() {
            if (this.cancelled) {
                return;
            }
            started = true;

            try {
                this.result = doDownload();
            } catch (IOException e) {
                log.error(e.getMessage(), e);
            } finally {
                this.done();
            }

        }

        /**
         * Return the result. Must be called after run is complete
         *
         * @return
         */
        public RunnableResult getResult() {
            if (!this.done) throw new IllegalStateException("Must wait for run to finish before getting result");
            return this.result;
        }

        private RunnableResult doDownload() throws IOException {

            log.info("Downloading " + srcUrl + " to " + outputFile.getAbsolutePath());

            HttpURLConnection conn = openConnection(this.srcUrl, null);

            long contentLength = -1;
            String contentLengthString = conn.getHeaderField("Content-Length");
            if (contentLengthString != null) {
                contentLength = Long.parseLong(contentLengthString);
            }

            InputStream is = null;
            OutputStream out = null;

            long downloaded = 0;
            long downSinceLast = 0;
            String curStatus;
            String msg1 = String.format("downloaded of %s total", contentLength >= 0 ? bytesToByteCountString(contentLength) : "unknown");
            int perc = 0;
            try {
                is = getInputStream(conn);
                out = new FileOutputStream(outputFile);

                byte[] buf = new byte[64 * 1024];
                int counter = 0;
                int interval = 100;
                int bytesRead = 0;
                while (!this.cancelled && (bytesRead = is.read(buf)) != -1) {
                    out.write(buf, 0, bytesRead);
                    downloaded += bytesRead;
                    downSinceLast += bytesRead;
                    counter = (counter + 1) % interval;
                    if (counter == 0 && this.monitor != null) {
                        curStatus = String.format("%s %s", bytesToByteCountString(downloaded), msg1);
                        this.monitor.updateStatus(curStatus);
                        if (contentLength >= 0) {
                            perc = (int) ((downSinceLast * 100) / contentLength);
                            this.monitor.fireProgressChange(perc);
                            if (perc >= 1) downSinceLast = 0;
                        }
                    }
                }
                log.info("Download complete.  Total bytes downloaded = " + downloaded);
            } catch (IOException e) {
                readErrorStream(conn);
                throw e;
            } finally {
                if (is != null) is.close();
                if (out != null) {
                    out.flush();
                    out.close();
                }
            }
            long fileLength = outputFile.length();

            if (this.cancelled) return RunnableResult.CANCELLED;

            boolean knownComplete = contentLength == fileLength;
            //Assume success if file length not known
            if (knownComplete || contentLength < 0) {
                if (this.monitor != null) {
                    this.monitor.fireProgressChange(100);
                    this.monitor.updateStatus("Done");
                }
                return RunnableResult.SUCCESS;
            } else {
                return RunnableResult.FAILURE;
            }

        }

        protected void done() {
            this.done = true;
        }

        public boolean isDone() {
            return this.done;
        }

        /**
         * See {@link java.util.concurrent.FutureTask#cancel(boolean)}
         *
         * @param mayInterruptIfRunning
         * @return
         */
        public boolean cancel(boolean mayInterruptIfRunning) {
            if (this.started && !mayInterruptIfRunning) {
                return false;
            }
            this.cancelled = true;
            return true;
        }

        public void setMonitor(ProgressMonitor monitor) {
            this.monitor = monitor;
        }

        /**
         * Convert bytes to human-readable string.
         * e.g. 102894 -> 102.89 KB. If too big or too small,
         * doesn't append a prefix just returns {@code bytes + " B"}
         *
         * @param bytes
         * @return
         */
        public String bytesToByteCountString(long bytes) {
            int unit = 1000;
            String prefs = "KMGT";

            if (bytes < unit) return bytes + " B";
            int exp = (int) (Math.log(bytes) / Math.log(unit));
            if (exp <= 0 || exp >= prefs.length()) return bytes + " B";

            String pre = (prefs).charAt(exp - 1) + "";
            return String.format("%.2f %sB", bytes / Math.pow(unit, exp), pre);
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


    /**
     * Extension of CookieManager that grabs cookies from the GenomeSpace identity server to store locally.
     * This is to support the GenomeSpace "single sign-on". Examples ...
     * gs-username=igvtest; Domain=.genomespace.org; Expires=Mon, 21-Jul-2031 03:27:23 GMT; Path=/
     * gs-token=HnR9rBShNO4dTXk8cKXVJT98Oe0jWVY+; Domain=.genomespace.org; Expires=Mon, 21-Jul-2031 03:27:23 GMT; Path=/
     */

    static class IGVCookieManager extends CookieManager {

        @Override
        public Map<String, List<String>> get(URI uri, Map<String, List<String>> requestHeaders) throws IOException {
            Map<String, List<String>> headers = new HashMap<String, List<String>>();
            headers.putAll(super.get(uri, requestHeaders));


            if (GSUtils.isGenomeSpace(uri.toURL())) {
                String token = GSUtils.getGSToken();
                if (token != null) {
                    List<String> cookieList = headers.get("Cookie");
                    boolean needsTokenCookie = true;
                    if (cookieList == null) {
                        cookieList = new ArrayList<String>(1);
                        headers.put("Cookie", cookieList);
                    }

                    for (String cookie : cookieList) {
                        if (cookie.startsWith("gs-token")) {
                            needsTokenCookie = false;
                            break;
                        }
                    }
                    if (needsTokenCookie) {
                        cookieList.add("gs-token=" + token);
                    }
                }
            }

            return Collections.unmodifiableMap(headers);
        }

        @Override
        public void put(URI uri, Map<String, List<String>> responseHeaders) throws IOException {
            String urilc = uri.toString().toLowerCase();
            if (urilc.contains("identity") && urilc.contains("genomespace")) {
                List<String> cookies = responseHeaders.get("Set-Cookie");
                if (cookies != null) {
                    for (String cstring : cookies) {
                        List<HttpCookie> cookieList = HttpCookie.parse(cstring);
                        for (HttpCookie cookie : cookieList) {
                            String cookieName = cookie.getName();
                            String value = cookie.getValue();
                            if (cookieName.equals("gs-token")) {
                                //log.debug("gs-token: " + value);
                                GSUtils.setGSToken(value);
                            } else if (cookieName.equals("gs-username")) {
                                //log.debug("gs-username: " + value);
                                GSUtils.setGSUser(value);
                            }
                        }
                    }
                }
            }
            super.put(uri, responseHeaders);
        }
    }

}
