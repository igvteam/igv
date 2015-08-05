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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv;

import org.apache.log4j.Logger;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.Hashtable;
import java.util.Map;

/**
 * An alternative to the preference API.  Introduced to work around limits on
 * string length in the default Linux implementation.
 * of Preferences
 *
 * @author jrobinso
 */
public class IGVPreferences {

    private static Logger log = Logger.getLogger(IGVPreferences.class);

    /**
     * Cache for preferences which should persist for this session only
     */
    static Hashtable<String, String> sessionCache = new Hashtable();

    static Hashtable<String, String> aliasTable = new Hashtable<String, String>();

    static {
        aliasTable.put("SAM>SORT_OPTION", "SAM.SORT_OPTION");
        aliasTable.put("FLAKING_REGIONS", "FLANKING_REGION");
    }

    Hashtable<String, String> userPreferences = null;

    private File prefFile;

    public IGVPreferences() {
        this.prefFile = null;
    }

    public IGVPreferences(File prefFile) {
        this.prefFile = prefFile;
    }

    public void put(String key, String value) {

        // Remove from session only, explicitly setting this overrides
        sessionCache.remove(key);

        if (aliasTable.containsKey(key)) {
            key = aliasTable.get(key);
        }
        if (userPreferences == null) {
            loadUserPreferences();
        }
        if (value == null) {
            userPreferences.remove(key);
        } else {
            userPreferences.put(key, value);
        }
        storePreferences();
    }

    public void putOverride(String key, String value) {
        sessionCache.put(key, value);
    }

    public String get(String key, String defaultValue) {
        if (sessionCache.containsKey(key)) {
            return sessionCache.get(key);
        }
        String val = get(key);
        return val == null ? defaultValue : val;
    }


    private String get(String key) {
        if (userPreferences == null) {
            loadUserPreferences();
        }
        return userPreferences.get(key);
    }

    public void remove(String key) {
        sessionCache.remove(key);
        userPreferences.remove(key);
        storePreferences();
    }

    public void clear() {
        userPreferences.clear();
        storePreferences();

    }

    synchronized void loadUserPreferences() {
        userPreferences = new Hashtable();
        if (prefFile == null) {
            prefFile = DirectoryManager.getPreferencesFile();
        }
        if (prefFile.exists()) {
            String prefFileName = prefFile.getAbsolutePath();
            load(prefFileName, false);
        }
    }

    /**
     * Load an override preference file.  These preferences to not persist as user preferences.
     *
     * @param path
     */
    public void loadOverrides(String path) {
        if (userPreferences == null) {
            loadUserPreferences();
        }
        load(path, true);
    }

    /**
     * @param prefFileName either a comma delimited list of key-value pairs, or a path to a properties file
     * @param override
     */
    private void load(String prefFileName, boolean override) {

        if (prefFileName.contains("=")) {
            String[] kvPairs = prefFileName.split(",");
            for (String kvPair : kvPairs) {
                String[] kv = kvPair.split("=");
                if (kv.length == 2) {
                    String key = kv[0];
                    String value = kv[1];
                    override(key, value, override);
                }
            }

        } else {
            BufferedReader reader = null;
            try {
                reader = ParsingUtils.openBufferedReader(prefFileName);
                String nextLine = null;
                while ((nextLine = reader.readLine()) != null) {
                    int idx = nextLine.indexOf('=');
                    if (idx > 0) {
                        String key = nextLine.substring(0, idx);
                        String value = nextLine.substring(idx + 1);
                        override(key, value, override);
                    }
                }
            } catch (IOException e) {
                log.error("Error loading preferences", e);
            } finally {
                try {
                    if (reader != null) {
                        reader.close();
                    }
                } catch (IOException ex) {
                    log.error("Error closing preferences file", ex);
                }
            }
        }
    }

    private void override(String key, String value, boolean override) {
        if (aliasTable.containsKey(key)) {
            key = aliasTable.get(key);
        }
        if (!value.equals("null")) {
            if (override) {
                log.info("Overriding preference: " + key + "=" + value);
                sessionCache.put(key, value);
            } else {
                if (userPreferences == null) {
                    loadUserPreferences();
                }
                userPreferences.put(key, value);
            }
        }
    }


    private synchronized void storePreferences() {

        if (userPreferences != null) {
            FileWriter fileWriter = null;
            try {
                fileWriter = new FileWriter(prefFile);
                PrintWriter pw = new PrintWriter(new BufferedWriter(fileWriter));
                for (Map.Entry<String, String> entry : userPreferences.entrySet()) {
                    pw.print(entry.getKey());
                    pw.print("=");
                    pw.println(entry.getValue());
                }
                pw.flush();
                pw.close();
            } catch (IOException e) {
                log.error("Error loading preferences", e);
            } finally {

                if (fileWriter != null) {
                    try {
                        fileWriter.close();
                    } catch (IOException e) {
                        // Ignore
                    }
                }
            }

        }
    }
}
