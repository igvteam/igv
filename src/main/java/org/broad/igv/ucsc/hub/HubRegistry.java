package org.broad.igv.ucsc.hub;

import com.google.gson.Gson;
import org.broad.igv.DirectoryManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.net.URL;
import java.util.*;

public class HubRegistry {

    public static final String UCSC_REST_PUBLICHUBS = "https://api.genome.ucsc.edu/list/publicHubs";
    public static final String BACKUP_HUBS_URL = "https://raw.githubusercontent.com/igvteam/igv-genomes/refs/heads/main/hubs/ucsc/publicHubs.json";
    public static final String UCSC_GENOMES = "https://api.genome.ucsc.edu/list/ucscGenomes";
    public static final String UCSC_GENOMES_BACKUP = "https://raw.githubusercontent.com/igvteam/igv-genomes/refs/heads/main/hubs/ucsc/ucscGenomes.json";


    private static Logger log = LogManager.getLogger(HubRegistry.class);


    private static List<HubDescriptor> allHubDescriptors = null;
    private static List<HubDescriptor> selectedHubDescriptors = null;
    private static Map<String, List<HubDescriptor>> selectedHubsMap = null;
    private static Set<String> ucscGenomeIDs = null;

    private HubRegistry() {
        // Private constructor to prevent instantiation
    }

    /**
     * Get all available hubs. This includes hubs loaded from  from the UCSC REST API and directly by URL.
     *
     * @return List of HubDescriptor objects
     */
    public static List<HubDescriptor> getAllHubs() {
        if (allHubDescriptors == null) {
            allHubDescriptors = loadDescriptors(UCSC_REST_PUBLICHUBS);
        }
        return allHubDescriptors;
    }

    /**
     * Get all user selected hubs. This is a subset of all hubs.
     *
     * @return List of HubDescriptor objects
     */
    public static List<HubDescriptor> getAllSelectedHubs() {
        return selectedHubDescriptors;
    }

    /**
     * Get all user selected hubs for a specific genome.
     *
     * @param ucscId
     * @return List of HubDescriptor objects
     */
    public static List<HubDescriptor> getSelectedHubsForGenome(String ucscId) {
        if (selectedHubsMap == null) {
            initializeSelectedHubsMap();
        }
        return selectedHubsMap.get(ucscId);
    }

    public static void addUserHub(Hub hub) {
        if (selectedHubDescriptors == null) {
            selectedHubDescriptors = new ArrayList<>();
        }
        selectedHubDescriptors.add(hub.getDescriptor());
        resetSelectedHubsMap();
        writeSelectedHubs();
    }

    /**
     * Remove an individual hub.
     *
     * @param hubDescriptor
     */
    public static void removeHub(HubDescriptor hubDescriptor) {
        selectedHubDescriptors.removeIf(hd -> hd.getUrl().equals(hubDescriptor.getUrl()));
        resetSelectedHubsMap();
        writeSelectedHubs();
    }

    /**
     * Replace the list of selected hubs with the provided list. This is called from the selection dialog
     *
     * @param hubDescriptors
     */
    public static void setSelectedHubs(List<HubDescriptor> hubDescriptors) {
        selectedHubDescriptors = hubDescriptors;
        resetSelectedHubsMap();
        writeSelectedHubs();
    }

    private static void writeSelectedHubs() {
        File file = new File(DirectoryManager.getIgvDirectory(), "hubs.txt");
        try (PrintWriter writer = new PrintWriter(file)) {
            // Write header
            writer.println("hubUrl\tshortLabel\tlongLabel\tdbList\tdescriptionUrl");

            // Write each hub descriptor
            for (HubDescriptor hub : selectedHubDescriptors) {
                writer.printf("%s\t%s\t%s\t%s\t%s%n",
                        hub.getUrl(),
                        hub.getShortLabel(),
                        hub.getLongLabel(),
                        hub.getDbList(),
                        hub.getDescriptionUrl());
            }
        } catch (IOException e) {
            log.error("Error writing selected hubs to file", e);
        }
    }

    private static List<HubDescriptor> readSelectedHubs() {
        selectedHubDescriptors = new ArrayList<>();
        File file = new File(DirectoryManager.getIgvDirectory(), "hubs.txt");

        if (file.exists()) {
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)))) {
                // Skip the header line
                String line = reader.readLine();

                // Read each line and create HubDescriptor objects
                while ((line = reader.readLine()) != null) {
                    String[] fields = line.split("\t");
                    if (fields.length == 5) {
                        String url = fields[0];
                        String shortLabel = fields[1];
                        String longLabel = fields[2];
                        String dbList = fields[3];
                        String descriptionUrl = fields[4];
                        selectedHubDescriptors.add(new HubDescriptor(url, shortLabel, longLabel, dbList, descriptionUrl));
                    }
                }
            } catch (IOException e) {
                log.error("Error reading selected hubs from file", e);
            }
        } else {
            log.warn("Selected hubs file does not exist: " + file.getAbsolutePath());
        }

        return selectedHubDescriptors;
    }

    private static void initializeSelectedHubsMap() {
        selectedHubsMap = new HashMap<>();
        selectedHubDescriptors = readSelectedHubs();
        resetSelectedHubsMap();
    }

    private static void resetSelectedHubsMap() {
        selectedHubsMap.clear();
        for (HubDescriptor hub : selectedHubDescriptors) {
            String dbList = hub.getDbList();
            String[] dbs = dbList.split(",");
            for (String db : dbs) {
                db = db.trim();
                if (selectedHubsMap.containsKey(db)) {
                    selectedHubsMap.get(db).add(hub);
                } else {
                    List<HubDescriptor> hubList = new ArrayList<>();
                    hubList.add(hub);
                    selectedHubsMap.put(db, hubList);
                }
            }
        }
    }

    public static List<HubDescriptor> loadDescriptors(String apiUrl) {

        List<HubDescriptor> hubDescriptors = new ArrayList<>();

        Map<String, Object> jsonResponse = null;
        Gson gson = new Gson();
        try {
            URL url = new URL(apiUrl);

            String jsonString = HttpUtils.getInstance().getContentsAsJSON(url);
            jsonResponse = gson.fromJson(jsonString, Map.class);

        } catch (Exception e) {
            // Try backup URL
            try {
                log.error("Failed to load Hub JSON.  Trying backup URL", e);
                String jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(BACKUP_HUBS_URL));
                jsonResponse = gson.fromJson(jsonString, Map.class);
            } catch (IOException ex) {
                log.error("Failed to load Hub JSON", e);
            }
        }

        if (jsonResponse != null) {
            Object hubObj = jsonResponse.get("publicHubs");
            for (Map<String, Object> hub : (List<Map<String, Object>>) hubObj) {
                String hubUrl = hub.get("hubUrl").toString();
                String shortLabel = hub.get("shortLabel").toString();
                String longLabel = hub.get("longLabel").toString();
                String dbList = hub.get("dbList").toString();
                String descriptionUrl = hub.get("descriptionUrl").toString();
                hubDescriptors.add(new HubDescriptor(hubUrl, shortLabel, longLabel, dbList, descriptionUrl));
            }
        }

        return hubDescriptors;
    }


    private Set<String> getUcscGenomeIDs() {
        if (ucscGenomeIDs == null) {
            try {
                String jsonString = getContentsAsJSON(UCSC_GENOMES, UCSC_GENOMES_BACKUP);
                Gson gson = new Gson();
                Map<String, Object> jsonResponse = gson.fromJson(jsonString, Map.class);
                Map<String, Object> hubGenomes = (Map<String, Object>) jsonResponse.get("ucscGenomes");
                ucscGenomeIDs = hubGenomes.keySet();
            } catch (IOException e) {
                log.error("Error loading UCSC hub genomes", e);
            }
        }
        return ucscGenomeIDs;
    }

    private String getContentsAsJSON(String url, String backupURL) throws IOException {
        String jsonString = null;
        try {
            jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(url));
        } catch (IOException e) {
            log.error("Failed to load JSON from " + url + ". Trying backup URL");
            jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(backupURL));
        }
        return jsonString;
    }

}

