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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.*;
import java.net.URL;

/**
 * 
 * This is a mostly deprecated class, created before range-byte requests were implemented.  It is still neccessary for a few
 * organizations, chiefly the Partners network, that inexplicably strip range-byte headers off outgoing http requests.
 * 
 * @author jrobinso
 */
public class SequenceServletWrapper {

    /**
     * Field description
     */
    public static boolean SEQUENCE_SERVER_AVAILALBLE = true;
    public static final int CONNECTION_TIMEOUT = 20000;
    private static Logger logger = Logger.getLogger(SequenceServletWrapper.class);


    public static byte[] readBytes(String urlString, String chr, int start, int end) {

        byte[] bytes = new byte[end - start];

        if (!SEQUENCE_SERVER_AVAILALBLE) {
            return bytes;
        }

        InputStream cis = null;
        try {

            URL url = new URL(urlString + "?chr=" + chr + "&start=" + start + "&end=" + end);

            cis = HttpUtils.getInstance().openConnectionStream(url);

            DataInputStream is = new DataInputStream(new BufferedInputStream(cis));
            int offset = 0;
            int numRead;
            while ((offset < bytes.length) && (numRead = is.read(bytes, offset, bytes.length - offset)) >= 0) {
                offset += numRead;
            }
            is.close();

            return bytes;


        } catch (IOException ex) {

            // Log connection errors once per session to prevent filling the log
            if (SEQUENCE_SERVER_AVAILALBLE) {
                SEQUENCE_SERVER_AVAILALBLE = false;

                showUnavailableMessage();

                logger.error("Error retrieving sequence from : " + urlString + ex.getMessage());
            }
            return null;
        }
        finally {
            if(cis != null) try {
                cis.close();
            } catch (IOException e) {
                
            }
        }
    }

    private static void showUnavailableMessage() throws HeadlessException {
        MessageUtils.showMessage(
                "<html>The IGV server at the Broad Institute is currently unavailable.  " +
                        "Features that require a reference sequence, <br>" +
                        "such as displaying alignment mismatches, will be disabled.  " +
                        "If this problem persists please <br>" +
                        "send email to igv-help@broadinstitute.org");
    }


}
