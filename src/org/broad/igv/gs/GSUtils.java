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

package org.broad.igv.gs;


import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 * @date Jun 9, 2011
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
        if (gsToken == null || !gsToken.equals(newToken)) {
            gsToken = newToken;

            if (Globals.isTesting()) {
                return;
            }

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
        if (gsToken == null && !Globals.isTesting()) {
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
        if (gsUser == null || !gsUser.equals(newUser)) {
            gsUser = newUser;

            if (Globals.isTesting()) {
                return;
            }

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
        if (gsUser == null && !Globals.isTesting()) {
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
        //String dmServer = getDMHost();
        //String atmServer = getATMHost();
        return url.getHost().contains("genomespace.org");
    }

    private synchronized static String getDMHost() {
        try {
            return (new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER))).getHost();
        } catch (MalformedURLException e) {
            log.error(e);
            return null;
        }
    }


    private synchronized static String getATMHost() {

        try {
            return (new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER))).getHost();
        } catch (MalformedURLException e) {
            log.error(e);
            return null;
        }

    }

    /**
     * Convert GenomeSpace format string to igv format type (extension)
     * Example: https://dmtest.genomespace.org:8444/datamanager/files/users/SAGDemo/Step1/TF.data.tab
     * ?dataformat=http://www.genomespace.org/datamanager/dataformat/gct/0.0.0
     *
     * @param format
     * @return
     */
    public static String parseDataFormatString(String format) {
        int versionIdx = format.lastIndexOf("/");
        String tmp = format.substring(0, versionIdx);
        int extIdx = tmp.lastIndexOf("/") + 1;
        return "." + tmp.substring(extIdx);
    }
}