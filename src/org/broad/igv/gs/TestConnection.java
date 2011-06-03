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

package org.broad.igv.gs;

import org.apache.commons.httpclient.Credentials;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.UsernamePasswordCredentials;
import org.apache.commons.httpclient.auth.AuthScope;
import org.apache.commons.httpclient.methods.GetMethod;
import org.broad.igv.util.IGVHttpUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;

/**
 * @author jrobinso
 * @date Jun 3, 2011
 */
public class TestConnection {

    

    public static void main(String [] args) throws IOException {

        URL url = new URL("https://dmtest.genomespace.org:8444/datamanager/files/users/ted/SubDirExample/all_aml_train_filt.gct");
        HttpClient client = new HttpClient();
        Credentials defaultcreds = new UsernamePasswordCredentials("ted", "qw");
        client.getState().setCredentials(new AuthScope("identitytest.genomespace.org", 8443, AuthScope.ANY_REALM), defaultcreds);

        GetMethod get = new GetMethod(url.toExternalForm());
        get.setDoAuthentication(true);

        client.getParams().setParameter("http.protocol.allow-circular-redirects", true);

        int status = client.executeMethod(get);
        System.out.println("status = " + status);

        InputStream is = get.getResponseBodyAsStream();

        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String nextLine;
        while((nextLine = br.readLine()) != null) {
            System.out.println(nextLine);
        }

    }

}
