package org.broad.igv.ucsc.hub;

import org.broad.igv.DirectoryManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.stream.Collectors;

public class HubRegistry {

    public static final String UCSC_REST_PUBLICHUBS = "https://api.genome.ucsc.edu/list/publicHubs";
    public static final String BACKUP_HUBS_URL = "https://raw.githubusercontent.com/igvteam/igv-genomes/refs/heads/main/hubs/ucsc/publicHubs.json";


    private static Logger log = LogManager.getLogger(HubRegistry.class);


    private static List<HubDescriptor> selectedHubDescriptors = null;
    private static List<HubDescriptor> ucscPublicHubDescriptors = null;

    private static Map<String, List<HubDescriptor>> allHubsMap = null;
    private static Map<String, List<HubDescriptor>> selectedHubsMap = null;

    private HubRegistry() {
        // Private constructor to prevent instantiation
    }

    /**
     * Get all user selected hubs.
     *
     * @return List of HubDescriptor objects
     */
    private static List<HubDescriptor> getSelectedHubs() {
        if (selectedHubDescriptors == null) {
            selectedHubDescriptors = readSelectedHubs();
        }
        return selectedHubDescriptors;
    }

    private static List<HubDescriptor> getUcscPublicHubDescriptors() {
        if (ucscPublicHubDescriptors == null) {
            ucscPublicHubDescriptors = loadUCSCDescriptors();
        }
        return ucscPublicHubDescriptors;
    }

    /**
     * Get all hubs (UCSC public + user selected) for a specific genome.
     *
     * @param ucscID  -- The UCSC genome ID, which may be different from the IGV genome ID.   For example hg38-1kg => hg38
     * @return List of HubDescriptor objects
     */
    public static List<HubDescriptor> getAllHubsForGenome(String ucscID) {
        if (allHubsMap == null) {
            List<HubDescriptor> allHubDescriptors = new ArrayList<>(getUcscPublicHubDescriptors());

            // Add user-selected hubs not already in the list
            Set<String> existingUrls = allHubDescriptors.stream()
                    .map(HubDescriptor::getUrl)
                    .collect(Collectors.toSet());
            List<HubDescriptor> hubsToAdd = getSelectedHubs().stream()
                    .filter(hub -> !existingUrls.contains(hub.getUrl()))
                    .collect(Collectors.toList());
            allHubDescriptors.addAll(0, hubsToAdd);

            // Build the map of hubs by genome
            allHubsMap = allHubDescriptors.stream()
                    .flatMap(hub -> Arrays.stream(hub.getDbList().split(","))
                            .map(db -> Map.entry(db.trim(), hub)))
                    .collect(Collectors.groupingBy(
                            Map.Entry::getKey,
                            Collectors.mapping(Map.Entry::getValue, Collectors.toList())
                    ));
        }

        return allHubsMap.getOrDefault(ucscID, Collections.emptyList());
    }


    /**
     * Get all user selected hubs for a specific genome.
     *
     * @param ucscID  -- The UCSC genome ID, which may be different from the IGV genome ID.   For example hg38-1kg => hg38
     * @return List of HubDescriptor objects
     */
    public static List<HubDescriptor> getSelectedHubsForGenome(String ucscID) {
        if (selectedHubsMap == null) {
            selectedHubsMap = getSelectedHubs().stream()
                    .flatMap(hub -> Arrays.stream(hub.getDbList().split(","))
                            .map(db -> Map.entry(db.trim(), hub)))
                    .collect(Collectors.groupingBy(
                            Map.Entry::getKey,
                            Collectors.mapping(Map.Entry::getValue, Collectors.toList())
                    ));
        }
        return selectedHubsMap.getOrDefault(ucscID, Collections.emptyList());
    }

    public static void addUserHub(Hub hub) {
        if (selectedHubDescriptors == null) {
            selectedHubDescriptors = new ArrayList<>();
        }
        selectedHubDescriptors.add(hub.getDescriptor());
        writeSelectedHubs();

        // Reset the maps so they will be rebuilt on next access
        selectedHubsMap = null;
        allHubsMap = null;
    }


    /**
     * Replace the list of selected hubs with the provided list. This is called from the selection dialog
     *
     * @param hubDescriptors
     */
    public static void setSelectedHubs(List<HubDescriptor> hubDescriptors) {
        selectedHubDescriptors = new ArrayList<>(hubDescriptors);
        writeSelectedHubs();

        // Reset the maps so they will be rebuilt on next access
        selectedHubsMap = null;
        allHubsMap = null;
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

    /**
     * Read the list of user selected hubs from a file.  These can include hubs selected from the UCSC list
     * as well as hubs added directly by URL.
     *
     * @return List of HubDescriptor objects
     */
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


    /**
     * Load UCSC hub descriptors from the UCSC API URL.
     *
     * @return List of HubDescriptor objects
     */
    public static List<HubDescriptor> loadUCSCDescriptors() {
        List<HubDescriptor> hubDescriptors = new ArrayList<>();
        org.json.JSONObject jsonResponse = fetchJsonResponse(UCSC_REST_PUBLICHUBS, BACKUP_HUBS_URL);

        if (jsonResponse != null) {
            org.json.JSONArray hubArray = jsonResponse.optJSONArray("publicHubs");
            if (hubArray != null) {
                for (int i = 0; i < hubArray.length(); i++) {
                    org.json.JSONObject hub = hubArray.getJSONObject(i);
                    hubDescriptors.add(new HubDescriptor(
                            hub.optString("hubUrl"),
                            hub.optString("shortLabel"),
                            hub.optString("longLabel"),
                            hub.optString("dbList"),
                            hub.optString("descriptionUrl")
                    ));
                }
            }
        }
        return hubDescriptors;
    }

    private static org.json.JSONObject fetchJsonResponse(String primaryUrl, String backupUrl) {
        try {
            String jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(primaryUrl));
            return new org.json.JSONObject(jsonString);
        } catch (Exception e) {
            log.error("Failed to load public hub list from " + primaryUrl, e);
            try {
                String jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(backupUrl));
                return new org.json.JSONObject(jsonString);
            } catch (IOException ex) {
                log.error("Failed to load public hub list from " + backupUrl, ex);
            }
        }
        return null;
    }

}

