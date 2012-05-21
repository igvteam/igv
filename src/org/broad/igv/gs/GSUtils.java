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

package org.broad.igv.gs;


import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * @author jrobinso
 * @date Jun 9, 2011
 *
 */
public class GSUtils {
    static final Logger log = Logger.getLogger(GSUtils.class);


    /*
    * Directory and filenames to save the token and username to facilitate SSO
    */
    private static String tokenSaveDir = ".gs";
    private static String tokenSaveFileName = ".gstoken";
    private static String usernameSaveFileName = ".gsusername";
    public static String gsUser = null;
    public static String gsToken = null;

    private static File getTokenSaveDir() {
        String userDir = System.getProperty("user.home");
        File gsDir = new File(userDir, tokenSaveDir);
        if (!gsDir.exists()) {
            gsDir.mkdir();
        }
        return gsDir;
    }

    private static File getTokenFile() {
        File gsDir = getTokenSaveDir();
        return (gsDir != null && gsDir.exists()) ? new File(gsDir, tokenSaveFileName) : null;
    }

    private static File getUsernameFile() {
        File gsDir = getTokenSaveDir();
        return (gsDir != null && gsDir.exists()) ? new File(gsDir, usernameSaveFileName) : null;
    }

    public static void setGSToken(String newToken) {
        if (Globals.isTesting()) {
            return;
        }
        if (gsToken == null || !gsToken.equals(newToken)) {
            gsToken = newToken;
            BufferedWriter bw = null;

            File gsDir = getTokenSaveDir();
            if (!gsDir.isDirectory()) {
                log.error("Could not store token for SSO.  File " + gsDir.getAbsolutePath() + "exists and is not a directory.");
                return; // someone made a file with this name...
            }
            File tokenFile = getTokenFile();
            if (tokenFile.exists()) tokenFile.delete();
            writeToFile(gsToken, tokenFile);
        }
    }

    public static String getGSToken() {
        if (Globals.isTesting()) {
            return null;
        }
        if (gsToken == null) {
            File file = GSUtils.getTokenFile();
            if (file.exists()) {
                BufferedReader br = null;
                try {
                    br = new BufferedReader(new FileReader(file));
                    gsToken = br.readLine();
                } catch (IOException e) {
                    log.error("Error reading GS cookie", e);
                } finally {
                    if (br != null) try {
                        br.close();
                    } catch (IOException e) {
                        // Ignore
                    }
                }
            }
        }
        return gsToken;
    }


    public static void setGSUser(String newUser) {
        if (Globals.isTesting()) {
            return;
        }
        if (gsUser == null || !gsUser.equals(newUser)) {
            gsUser = newUser;

            File gsDir = getTokenSaveDir();
            if (!gsDir.isDirectory()) {
                log.error("Could not store token for SSO.  File " + gsDir.getAbsolutePath() + "exists and is not a directory.");
                return; // someone made a file with this name...
            }
            File userFile = getUsernameFile();
            if (userFile.exists()) userFile.delete();
            writeToFile(gsUser, userFile);
        }
    }


    public static String getGSUser() throws IOException {
        if (Globals.isTesting()) {
            return null;
        }
        if (gsUser == null) {
            BufferedReader br = null;
            try {
                File tokenFile = getUsernameFile();
                if (tokenFile.exists()) {
                    br = new BufferedReader(new FileReader(tokenFile));
                    gsUser = br.readLine().trim();
                }
            } finally {
                try {
                    if (br != null) br.close();
                } catch (Exception e) {
                }
            }
            return gsUser;
        }
        return gsUser;
    }


    public static void logout() {

        gsToken = null;
        gsUser = null;
        gsToken = null;
        File userfile = getUsernameFile();
        if (userfile.exists()) {
            userfile.delete();
        }
        File tokenFile = getTokenFile();
        if (tokenFile.exists()) {
            tokenFile.delete();
        }


        try {
            URI gsURI = new URI(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER));
            final CookieStore cookieStore = ((CookieManager) CookieManager.getDefault()).getCookieStore();
            List<HttpCookie> cookies = new ArrayList<HttpCookie>(cookieStore.get(gsURI));
            if (cookies != null) {
                for (HttpCookie cookie : cookies) {
                    final String name = cookie.getName();
                    if (name.equals("gs-token") || name.equals("gs-username")) {
                        cookieStore.remove(gsURI, cookie);
                    }
                }
            }
        } catch (URISyntaxException e) {
           log.error("Error creating GS URI", e);
        }
    }

    private static void writeToFile(String line, File aFile) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(aFile));
            bw.write(line);

            bw.close();
        } catch (Exception e) {
            log.error("Failed to save the token for later Single Sign on", e);
        } finally {
            try {
                if (bw != null) bw.close();
            } catch (Exception e) {
            }
        }
    }


    public static boolean isGenomeSpace(URL url) {
        return url.getHost().contains("genomespace");
    }

}