package org.broad.igv.ucsc.hub;

import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.gson.Gson;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Map;


public class HubRegistry {


    static Map<String, List<String>> hubURLMap = null;
    static Map<String, List<HubDescriptor>> hubDescriptorMap = null;

    public static boolean hasHubsForGenome(String ucscId) {
        Map<String, List<String>> hubURLMap = getHubURLs();
        List<String> hubUrls = hubURLMap.get(ucscId);
        return hubUrls != null && !hubUrls.isEmpty();
    }

    public static List<String> getHubUrlsForGenome(String ucscId) {
        Map<String, List<String>> hubURLMap = getHubURLs();
        return hubURLMap == null ? null : hubURLMap.get(ucscId);
    }

    private static Map<String, List<String>> getHubURLs() {

        if (hubURLMap == null) {

            String filePath = PreferencesManager.getPreferences().get(Constants.AUXILLARY_HUBS_URL);

            hubURLMap = new HashMap<>();
            try (BufferedReader br = ParsingUtils.openBufferedReader(filePath)) {
                String line;
                String currentGenomeId = null;
                List<String> currentURLList = null;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("#")) {
                        continue;
                    }
                    line = line.trim();
                    if (currentGenomeId == null) {
                        currentGenomeId = line;
                        currentURLList = new ArrayList<>();
                        hubURLMap.put(currentGenomeId, currentURLList);
                    } else if (line.length() == 0) {
                        currentGenomeId = null;
                    } else {
                        currentURLList.add(line);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return hubURLMap;
    }


    public static List<HubDescriptor> getPublicHubs(String ucscId) {
        if (hubDescriptorMap == null) {
            loadPublicDescriptors();
        }
        return hubDescriptorMap.get(ucscId);
    }

    public static void loadPublicDescriptors() {

        if(hubDescriptorMap == null) {
            hubDescriptorMap = new HashMap<>();
        }

        Map<String, Object> jsonResponse = loadJSON();

        if (jsonResponse != null) {
            Object hubObj = jsonResponse.get("publicHubs");

            for (Map<String, Object> hub : (List<Map<String, Object>>) hubObj) {
                String dbList = hub.get("dbList").toString();
                String hubUrl = hub.get("hubUrl").toString();
                String shortLabel = hub.get("shortLabel").toString();
                String longLabel = hub.get("longLabel").toString();

                HubDescriptor hubDescriptor = new HubDescriptor(shortLabel, longLabel, hubUrl);
                String[] dbs = dbList.split(",");
                for(String db : dbs) {
                    db = db.trim();
                    if (hubDescriptorMap.containsKey(db)) {
                        hubDescriptorMap.get(db).add(hubDescriptor);
                    } else {
                        List<HubDescriptor> hubList = new ArrayList<>();
                        hubList.add(hubDescriptor);
                        hubDescriptorMap.put(db, hubList);
                    }
                }
            }
        }
    }

    public static Map<String, Object> loadJSON() {
        String apiUrl = "https://api.genome.ucsc.edu/list/publicHubs";
        Gson gson = new Gson();
        try {
            URL url = new URL(apiUrl);
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("GET");
            connection.setRequestProperty("Accept", "application/json");

            if (connection.getResponseCode() == HttpURLConnection.HTTP_OK) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()));
                Map<String, Object> jsonResponse = gson.fromJson(reader, Map.class);
                reader.close();
                return jsonResponse;
            } else {
                System.err.println("Failed to fetch data from API. HTTP code: " + connection.getResponseCode());
                return null;
            }
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }

    public static void main(String[] args) throws IOException {

        List<HubDescriptor> descriptors = getPublicHubs("hg38");
        if (descriptors != null) {
            for (HubDescriptor descriptor : descriptors) {
                System.out.println("Hub URL: " + descriptor.getUrl());
               // System.out.println("Short Label: " + descriptor.getShortLabel());
              //   System.out.println("Long Label: " + descriptor.getLongLabel());
            }
        } else {
            System.out.println("No public hubs found for the specified genome.");
        }
    }
}

