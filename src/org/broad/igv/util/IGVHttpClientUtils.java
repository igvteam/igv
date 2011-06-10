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

import org.apache.http.HttpResponse;
import org.apache.http.auth.AuthScope;
import org.apache.http.auth.UsernamePasswordCredentials;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.BasicResponseHandler;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.log4j.Logger;
import org.broad.igv.gs.GSLoginDialog;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
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

    static DefaultHttpClient client = null;

    /**
     * Execute a get on the url and return the response stream.  It is the responsibility of the caller to
     * close the stream.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public static InputStream openGSConnectionStream(URL url) throws IOException {
        HttpResponse response = executeGet(url);
        return response.getEntity().getContent();
    }

    /**
     * Execute a get on the url and return the header field value.
     *
     * @param url
     * @return
     * @throws IOException
     */
    public static String getGSHeaderField(URL url, String key) throws IOException {
        HttpResponse response = executeGet(url);
        return response.getFirstHeader(key).getValue();
    }

    private static HttpResponse executeGet(URL url) throws IOException {

        if (client == null) client = new DefaultHttpClient();
        HttpGet getMethod = null;
        try {

            if (GSUtils.isGenomeSpace(url)) {
                GSUtils.checkForCookie(client, url);
            }

            getMethod = new HttpGet(url.toExternalForm());
            //getMethod.setDoAuthentication(true);
            client.getParams().setParameter("http.protocol.allow-circular-redirects", true);
            HttpResponse response = client.execute(getMethod);
            final int statusCode = response.getStatusLine().getStatusCode();
            if (statusCode == 401) {
                getMethod.abort();
                // Try again
                client.getCredentialsProvider().clear();

                login(url);


                return executeGet(url);
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
        String userpass = getGSUserPass(owner);
        if (userpass == null) {
            throw new RuntimeException("Access denied to: " + url.toString());
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

        GSLoginDialog dlg = new GSLoginDialog(owner);
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
