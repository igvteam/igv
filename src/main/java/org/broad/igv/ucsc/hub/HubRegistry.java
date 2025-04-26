package org.broad.igv.ucsc.hub;

import com.google.gson.Gson;
import org.broad.igv.DirectoryManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.*;

public class HubRegistry {

    public static final String UCSC_REST_PUBLICHUBS = "https://api.genome.ucsc.edu/list/publicHubs";

    private static Logger log = LogManager.getLogger(HubRegistry.class);


    static List<HubDescriptor> allHubDescriptors = null;

    static Map<String, List<HubDescriptor>> hubDescriptorMap = null;

    static Map<String, List<HubDescriptor>> selectedHubsMap = null;

    public static List<HubDescriptor> getAllHubDescriptors() {
        if (allHubDescriptors == null) {
            allHubDescriptors = loadDescriptors(UCSC_REST_PUBLICHUBS);
        }
        return allHubDescriptors;
    }

    public static void setSelectedHubs(List<HubDescriptor> hubDescriptors) {

        resetSelectedHubsMap(hubDescriptors);

        File file = new File(DirectoryManager.getIgvDirectory(), "hubs.txt");
        try (PrintWriter writer = new PrintWriter(file)) {
            // Write header
            writer.println("shortLabel\tlongLabel\turl\tdescriptionUrl\tdbList");

            // Write each hub descriptor
            for (HubDescriptor hub : hubDescriptors) {
                writer.printf("%s\t%s\t%s\t%s\t%s%n",
                        hub.getShortLabel(),
                        hub.getLongLabel(),
                        hub.getUrl(),
                        hub.getDescriptionUrl(),
                        hub.getDbList());
            }
        } catch (IOException e) {
            log.error("Error writing selected hubs to file", e);
        }
    }

    public static List<HubDescriptor> readSelectedHubs() {
        List<HubDescriptor> hubDescriptors = new ArrayList<>();
        File file = new File(DirectoryManager.getIgvDirectory(), "hubs.txt");

        if (file.exists()) {
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)))) {
                // Skip the header line
                String line = reader.readLine();

                // Read each line and create HubDescriptor objects
                while ((line = reader.readLine()) != null) {
                    String[] fields = line.split("\t");
                    if (fields.length == 5) {
                        String shortLabel = fields[0];
                        String longLabel = fields[1];
                        String url = fields[2];
                        String descriptionUrl = fields[3];
                        String dbList = fields[4];
                        hubDescriptors.add(new HubDescriptor(shortLabel, longLabel, url, descriptionUrl, dbList));
                    }
                }
            } catch (IOException e) {
                log.error("Error reading selected hubs from file", e);
            }
        } else {
            log.warn("Selected hubs file does not exist: " + file.getAbsolutePath());
        }

        return hubDescriptors;
    }

    public static List<HubDescriptor> getSelectedHubs(String ucscId) {
        if (selectedHubsMap == null) {
            initializeSelectedHubsMap();
        }
        return selectedHubsMap.get(ucscId);
    }

    public static List<HubDescriptor> getAllSelectedHubs() {
        if (selectedHubsMap == null) {
            initializeSelectedHubsMap();
        }
        List<HubDescriptor> allHubs = new ArrayList<>();
        for (List<HubDescriptor> hubs : selectedHubsMap.values()) {
            allHubs.addAll(hubs);
        }
        return allHubs;
    }

    private static void initializeSelectedHubsMap() {
        selectedHubsMap = new HashMap<>();
        List<HubDescriptor> hubDescriptors = readSelectedHubs();
        resetSelectedHubsMap(hubDescriptors);
    }

    private static void resetSelectedHubsMap(List<HubDescriptor> hubDescriptors) {
        selectedHubsMap.clear();
        for (HubDescriptor hub : hubDescriptors) {
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


    public static void addUserHub(Hub hub) {

    }


    public static List<HubDescriptor> loadDescriptors(String apiUrl) {

        List<HubDescriptor> hubDescriptors = new ArrayList<>();

        Map<String, Object> jsonResponse = null;
        Gson gson = new Gson();
        try {
            URL url = new URL(apiUrl);
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("GET");
            connection.setRequestProperty("Accept", "application/json");

            if (connection.getResponseCode() == HttpURLConnection.HTTP_OK) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()));
                jsonResponse = gson.fromJson(reader, Map.class);
                reader.close();

            } else {
                log.error("Failed to fetch HUB data from API " + apiUrl + ". HTTP code: " + connection.getResponseCode());

            }
        } catch (Exception e) {
            log.error("Failed to load Hub JSON", e);

        }

        if (jsonResponse != null) {
            Object hubObj = jsonResponse.get("publicHubs");
            for (Map<String, Object> hub : (List<Map<String, Object>>) hubObj) {
                String dbList = hub.get("dbList").toString();
                String hubUrl = hub.get("hubUrl").toString();
                String shortLabel = hub.get("shortLabel").toString();
                String longLabel = hub.get("longLabel").toString();
                String descriptionUrl = hub.get("descriptionUrl").toString();
                hubDescriptors.add(new HubDescriptor(shortLabel, longLabel, hubUrl, descriptionUrl, dbList));
            }
        }

        return hubDescriptors;
    }

    /**
     * Main method for testing purposes
     */
    public static void main(String[] args) throws IOException {

        List<HubDescriptor> descriptors = getAllHubDescriptors();
        if (descriptors != null) {
            for (HubDescriptor descriptor : descriptors) {
                //    System.out.println("Hub URL: " + descriptor.getUrl());
            }
        }
    }
}

