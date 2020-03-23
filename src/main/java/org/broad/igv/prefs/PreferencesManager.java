package org.broad.igv.prefs;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;

import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.util.ParsingUtils;


import java.io.*;
import java.util.*;
import java.util.List;

import static org.broad.igv.prefs.Constants.*;

/**
 * Created by jrobinso on 1/21/17.
 */
public class PreferencesManager implements IGVEventObserver {

    public static final String SEPARATOR_KEY = "---";
    public static final String INFO_KEY = "info";

    private static List<PreferenceGroup> preferenceGroupList;
    private static Logger log = Logger.getLogger(PreferencesManager.class);
    private static Map<String, IGVPreferences> preferencesMap = Collections.synchronizedMap(new HashMap<>());
    private static IGVPreferences genericDefaults;

    public static boolean forceDefaults = false;

    private static String prefFile;  // User preferences file

    static Hashtable<String, String> aliasTable = new Hashtable<String, String>();


    static {
        aliasTable.put("SAM>SORT_OPTION", "SAM.SORT_OPTION");
        aliasTable.put("FLAKING_REGIONS", "FLANKING_REGION");
    }


    private PreferencesManager() {

    }

    private static PreferencesManager theInstance = new PreferencesManager();

    public static IGVPreferences getPreferences(String category) {

        if (preferenceGroupList == null) {
            init();
        }
        if (forceDefaults) {
            return genericDefaults;
        } else if (preferencesMap.containsKey(category)) {
            return preferencesMap.get(category);
        } else {
            return preferencesMap.get(NULL_CATEGORY);
        }
    }

    private static void init() {

        try {
            preferenceGroupList = loadPreferenceList();

            Map<String, Map<String, String>> defaultPreferences = new HashMap<>();

            Map<String, String> nullCategory = loadDefaults23();

            defaultPreferences.put(NULL_CATEGORY, nullCategory);
            defaultPreferences.put(RNA, new HashMap<>());
            defaultPreferences.put(THIRD_GEN, new HashMap<>());

            for (PreferenceGroup group : preferenceGroupList) {
                String category = group.category == null ? NULL_CATEGORY : group.category;
                Map<String, String> defs = defaultPreferences.get(category);
                if (defs == null) {
                    log.info("Unknown preference category: " + category);
                }
                for (Preference pref : group.preferences) {
                    defs.put(pref.getKey(), pref.getDefaultValue());
                }
            }

            genericDefaults = new IGVPreferences(defaultPreferences.get(NULL_CATEGORY), null, null);

            Map<String, String> defaults = defaultPreferences.get(NULL_CATEGORY);
            Map<String, String> rnaDefaults = defaultPreferences.get(RNA);
            Map<String, String> thirdGenDefaults = defaultPreferences.get(THIRD_GEN);

            Map<String, Map<String, String>> userPrefs = loadUserPreferences();

            final IGVPreferences nullPrefs = new IGVPreferences(userPrefs.get(NULL_CATEGORY), defaults, null);
            extractMutationColors(nullPrefs);
            preferencesMap.put(NULL_CATEGORY, nullPrefs);
            preferencesMap.put(RNA, new IGVPreferences(userPrefs.get(RNA), rnaDefaults, nullPrefs));
            preferencesMap.put(THIRD_GEN, new IGVPreferences(userPrefs.get(THIRD_GEN), thirdGenDefaults, nullPrefs));


        } catch (IOException e) {
            e.printStackTrace();
        }

        IGVEventBus.getInstance().subscribe(PreferencesChangeEvent.class, theInstance);
    }

    private static void extractMutationColors(IGVPreferences prefs) {
        String cts = prefs.get("MUTATION_COLOR_TABLE");
        if(cts != null) {
            String [] tokens= cts.split(";");
            for(String t : tokens) {
                String [] kv = t.split("=");
                if(kv.length  == 2) {
                    String key = IGVPreferences.getMutationColorKey(kv[0]);
                    prefs.put(key, kv[1]);
                }
            }
            prefs.remove("MUTATION_COLOR_TABLE");
        }
    }

    public static IGVPreferences getPreferences() {
        return forceDefaults ? genericDefaults : getPreferences(NULL_CATEGORY);
    }

    public static Collection<IGVPreferences> getAllPreferences() {
        return preferencesMap.values();
    }

    public static void setPrefsFile(String prefsFile) {
        prefFile = prefsFile;
    }

    private static Map<String, Map<String, String>> loadUserPreferences() {

        try {
            if (prefFile == null) {
                prefFile = DirectoryManager.getPreferencesFile().getAbsolutePath();
            }
            return load(prefFile);
        } catch (Exception e) {
            log.error("Error loading preferences file: " + prefFile, e);
            return null;
        }

    }

    public static void updateAll(Map<String, Map<String, String>> preferenceMap) {
        for (Map.Entry<String, Map<String, String>> entry : preferenceMap.entrySet()) {
            IGVPreferences preferences = getPreferences(entry.getKey());
            if (preferences != null) {
                preferences.putAll(entry.getValue());
            }
        }
    }

    public static void loadOverrides(String overridePropertyFilePath) {

        if (preferencesMap.get(NULL_CATEGORY) == null) {
            loadUserPreferences();
        }

        Map<String, Map<String, String>> overrides = load(overridePropertyFilePath);
        for (Map.Entry<String, Map<String, String>> entry : overrides.entrySet()) {

            IGVPreferences prefs = preferencesMap.containsKey(entry.getKey()) ?
                    preferencesMap.get(entry.getKey()) :
                    preferencesMap.get(NULL_CATEGORY);

            prefs.addOverrides(entry.getValue());
        }
    }

    private static Map<String, Map<String, String>> load(String prefFileName) {

        Map<String, Map<String, String>> prefMap = new HashMap<>();
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(prefFileName);
            String nextLine = null;
            String category = NULL_CATEGORY;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("##")) {
                    category = nextLine.substring(2).trim();
                } else {
                    Map<String, String> prefs = prefMap.get(category);
                    if (prefs == null) {
                        prefs = Collections.synchronizedMap(new HashMap<>());
                        prefMap.put(category, prefs);
                    }
                    int idx = nextLine.indexOf('=');
                    if (idx > 0) {
                        KeyValue kv = translate(nextLine.substring(0, idx), nextLine.substring(idx + 1));
                        prefs.put(kv.key, kv.value);
                    }
                }
            }
        } catch (IOException e) {
            log.info("Error loading preferences " + e.getMessage());
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException ex) {
                log.error("Error closing preferences file", ex);
            }
        }
        return prefMap;
    }



    /**
     * Update legacy preference key/value
     * @param key
     * @param value
     * @return
     */
    private static KeyValue translate(String key, String value) {
        if (aliasTable.containsKey(key)) {
            key = aliasTable.get(key);
        }
        else if (key.equals(SAM_SHADE_BASES)) {
            boolean b = value.equalsIgnoreCase("quality") || value.equalsIgnoreCase("true");
            value = String.valueOf(b);
        }
        return new KeyValue(key, value);
    }



    public static List<PreferenceGroup> loadPreferenceList() throws IOException {

        List<PreferenceGroup> groupList = new ArrayList<>();
        List<Preference> prefList = null;

        BufferedReader reader = null;

        String nextLine;
        String group = null;

        try {
            reader = new BufferedReader(new InputStreamReader(PreferencesEditor.class.getResourceAsStream("/org/broad/igv/prefs/preferences.tab")));
            while ((nextLine = reader.readLine()) != null) {
                nextLine = nextLine.trim();

                if (nextLine.startsWith("//") || nextLine.length() == 0) {
                    continue;
                } else if (nextLine.startsWith(SEPARATOR_KEY)) {
                    prefList.add(new Preference(SEPARATOR_KEY, group));
                    continue;
                } else if (nextLine.startsWith(INFO_KEY)) {
                    Preference preference = new Preference(INFO_KEY, nextLine.substring(INFO_KEY.length()).trim(), group);
                    prefList.add(preference);   // "Blank" preference
                    continue;
                } else if (nextLine.startsWith("##")) {

                    group = null;  // End previous group
                    if (nextLine.length() > 2) {
                        group = nextLine.substring(2);  // New group
                    }
                    continue;
                } else if (nextLine.startsWith("#")) {

                    // New tab
                    String[] tokens = Globals.tabPattern.split(nextLine);
                    String tabLabel = tokens[0].substring(1);
                    String category = tokens.length > 1 ? tokens[1] : null;

                    prefList = new ArrayList<>();
                    PreferenceGroup preferenceGroup = new PreferenceGroup(tabLabel, category, prefList);
                    groupList.add(preferenceGroup);

                    group = null;

                    continue;
                } else {
                    String[] tokens = Globals.tabPattern.split(nextLine);
                    if (tokens.length < 3) {
                        if (tokens.length == 2) {
                            // Hidden preference (not shown in editor)
                            tokens = new String[]{tokens[0], "", "", tokens[1]};
                            prefList.add(new Preference(tokens, group));
                        }

                    } else {
                        prefList.add(new Preference(tokens, group));
                    }
                }
            }
        } finally {
            if (reader != null) reader.close();
        }

        return groupList;
    }

    /**
     * Load defaults form version 2.3.  Catch any preferences that were missed in the conversion.
     *
     * @return
     * @throws IOException
     */
    public static Map<String, String> loadDefaults23() throws IOException {

        Map<String, String> defs = new HashMap<>();

        BufferedReader reader = null;
        String nextLine;
        try {
            reader = new BufferedReader(new InputStreamReader(PreferencesEditor.class.getResourceAsStream("/org/broad/igv/prefs/defaults_2.3.tab")));
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = Globals.tabPattern.split(nextLine);
                if (tokens.length == 2) {
                    defs.put(tokens[0], tokens[1]);
                }
            }
        } finally {
            if (reader != null) reader.close();
        }
        return defs;
    }

    /**
     * Override a preference for this session.  We don't have a parameter to indicate experiment type so override
     * it for all preference categories.
     * @param prefKey
     * @param prefVal
     */
    public static void setOverride(String prefKey, String prefVal) {

        if (preferenceGroupList == null) {
            init();
        }
        for(IGVPreferences prefs : preferencesMap.values()) {
            prefs.override(prefKey, prefVal);
        }
    }

    private synchronized void storePreferences() {

        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(prefFile);
            PrintWriter pw = new PrintWriter(new BufferedWriter(fileWriter));

            for (Map.Entry<String, IGVPreferences> entry : preferencesMap.entrySet()) {
                if (!entry.getKey().equals(NULL_CATEGORY)) {
                    pw.println();
                    pw.println("##" + entry.getKey());
                }
                entry.getValue().print(pw);
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

    @Override
    public void receiveEvent(Object event) {
        if (event instanceof PreferencesChangeEvent) {
            storePreferences();
        }
    }

    static class KeyValue {
        String key;
        String value;

        public KeyValue(String key, String value) {
            this.key = key;
            this.value = value;
        }
    }

    static class Preference {

        String group;
        String[] tokens;

        Preference(String[] tokens, String group) {
            this.tokens = tokens;
            this.group = group;
        }

        Preference(String key, String group) {
            this(new String[]{key, null, null, null}, group);
        }

        Preference(String key, String label, String group) {
            this(new String[]{key, label, null, null}, group);
        }

        String getKey() {
            return tokens[0];
        }

        String getLabel() {
            return tokens[1];
        }

        String getType() {
            return tokens[2];
        }

        String getDefaultValue() {
            return tokens.length < 4 || tokens[3] == null || tokens[3].equals("null") ? null : tokens[3];
        }

        String getComment() {
            return tokens.length > 4 ? tokens[4] : null;
        }

        String getGroup() {
            return group;
        }

        String printString() {
            String str = getKey() + "\t" + getLabel() + "\t" + getType() + "\t" + getDefaultValue();
            if (getComment() != null) str += "\t" + getComment();
            return str;
        }
    }

    static class PreferenceGroup {

        String tabLabel;
        String category;
        List<Preference> preferences;

        public PreferenceGroup(String tabLabel, String category, List<Preference> preferences) {
            this.tabLabel = tabLabel;
            this.category = category;
            this.preferences = preferences;
        }
    }

}
