/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ga4gh;

import com.google.gson.*;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.*;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Helper class for  Google API prototype
 * <p/>
 * <p/>
 * Created by jrobinso on 8/15/14.
 */
public class Ga4ghAPIHelper {

    private static Logger log = Logger.getLogger(Ga4ghAPIHelper.class);

    public static final String RESOURCE_TYPE = "ga4gh";

    public static final Ga4ghProvider GA4GH_GOOGLE_PROVIDER = new Ga4ghProvider(
            "Google",
            "https://www.googleapis.com/genomics/v1beta",
            "AIzaSyC-dujgw4P1QvNd8i_c-I-S_P1uxVZzn0w",
            Arrays.asList(
                    new Ga4ghDataset("10473108253681171589", "1000 Genomes", "hg19"),
                    new Ga4ghDataset("383928317087", "PGP", "hg19"),
                    new Ga4ghDataset("461916304629", "Simons Foundation", "hg19")
            ));

    public static final Ga4ghProvider GA4GH_NCBI_PROVIDER = new Ga4ghProvider("NCBI", "http://trace.ncbi.nlm.nih.gov/Traces/gg", null,
            Arrays.asList(
                    new Ga4ghDataset("SRP034507", "SRP034507", "M74568"),
                    new Ga4ghDataset("SRP029392", "SRP029392", "NC_004917")
            ));

    public static final Ga4ghProvider[] providers = {
            //     new Ga4ghProvider("EBI", "http://193.62.52.16", null, Arrays.asList(new Ga4ghDataset("data", "data"))),
            GA4GH_GOOGLE_PROVIDER,
            GA4GH_NCBI_PROVIDER};


    final static Map<String, List<Ga4ghReadset>> readsetCache = new HashMap<String, List<Ga4ghReadset>>();


    public static List<Ga4ghReadset> readsetSearch(Ga4ghProvider provider, Ga4ghDataset dataset, int maxResults) throws IOException {

        String datasetId = dataset.getId();
        List<Ga4ghReadset> readsets = readsetCache.get(datasetId);

        if (readsets == null) {

            readsets = new ArrayList();

            String genomeId = genomeIdMap.get(provider.getName() + " " + datasetId); // Hack until meta data on readsets is available

            // Loop through pages
            int maxPages = 100;
            JsonPrimitive pageToken = null;
            while (maxPages-- > 0) {
                String contentToPost = "{" +
                        "\"datasetIds\": [\"" + datasetId + "\"]" +
                        (pageToken == null ? "" : ", \"pageToken\": " + pageToken) +
                        ", \"maxResults\":" + maxResults +
                        "}";

                String result = doPost(provider, "/readsets/search", contentToPost, null); //"fields=readsets(id,name, fileData),nextPageToken");

                if(result == null) return null;

                JsonParser parser = new JsonParser();
                JsonObject obj = parser.parse(result).getAsJsonObject();

                Iterator<JsonElement> iter = obj.getAsJsonArray("readsets").iterator();

                while (iter.hasNext()) {
                    JsonElement next = iter.next();
                    JsonObject jobj = next.getAsJsonObject();
                    String id = jobj.get("id").getAsString();
                    String name = jobj.get("name").getAsString();
                    readsets.add(new Ga4ghReadset(id, name, genomeId));
                }

                if (readsets.size() >= maxResults) break;

                pageToken = obj.getAsJsonPrimitive("nextPageToken");
                if (pageToken == null) break;
            }

            Collections.sort(readsets, new Comparator<Ga4ghReadset>() {
                @Override
                public int compare(Ga4ghReadset o1, Ga4ghReadset o2) {
                    return o1.getName().compareTo(o2.getName());
                }
            });

            readsetCache.put(datasetId, readsets);
        }

        return readsets;
    }


    public static List<Alignment> reads(Ga4ghProvider provider, String readsetId, String chr, int start, int end) throws IOException {

        List<Alignment> alignments = new ArrayList<Alignment>(10000);
        int maxPages = 10000;
        JsonPrimitive pageToken = null;
        StringBuffer result = new StringBuffer();

        while (maxPages-- > 0) {
            String contentToPost = "{" +
                    "\"readsetIds\": [\"" + readsetId + "\"]" +
                    ", \"sequenceName\": \"" + chr + "\"" +
                    ", \"sequenceStart\": \"" + start + "\"" +
                    ", \"sequenceEnd\": \"" + end + "\"" +
                    ", \"maxResults\": \"10000\"" +
                    // (pageToken == null ? "" : ", \"pageToken\": " + pageToken) +
                    "}";

            String readString = doPost(provider, "/reads/search", contentToPost, "");

            if (readString == null) {
                return null;
            }

            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(readString).getAsJsonObject();

            JsonArray reads = obj.getAsJsonArray("reads");

            Iterator<JsonElement> iter = reads.iterator();
            while (iter.hasNext()) {
                JsonElement next = iter.next();
                Ga4ghAlignment alignment = new Ga4ghAlignment(next.getAsJsonObject());
                alignments.add(alignment);
            }

            //System.out.println("# reads = " + reads.size());

            pageToken = obj.getAsJsonPrimitive("nextPageToken");
            if (pageToken == null) break;

        }

        //System.out.println("# pages= " + (10000 - maxPages));

        return alignments;

    }

    private static String doPost(Ga4ghProvider provider, String command, String content, String fields) throws IOException {

        String authKey = provider.getAuthKey();
        String baseURL = provider.getBaseURL();
        String token = OAuthUtils.getInstance().getAccessToken();

        String fullUrl = baseURL + command;
        if (authKey != null) {
            fullUrl += "?key=" + authKey;
        }
        if (fields != null) {
            fullUrl += (authKey == null ? "?" : "&") + fields;
        }

        URL url = new URL(fullUrl);


        byte[] bytes = content.getBytes();

        // Create a URLConnection
        HttpURLConnection connection = null;
        BufferedReader br;
        OutputStream outputStream;

        try {
            connection = (HttpURLConnection) url.openConnection();
            connection.setUseCaches(false);
            connection.setDoInput(true);
            connection.setDoOutput(true);
            connection.setRequestMethod("POST");
            //connection.setRequestProperty("Content-Length", "" + bytes.length);
            connection.setRequestProperty("Content-Type", "application/json");
            connection.setRequestProperty("Cache-Control", "no-cache");
            connection.setRequestProperty("Accept-Encoding", "gzip");
            connection.setRequestProperty("User-Agent", "IGV (gzip)");

            if (token != null) {
                connection.setRequestProperty("Authorization", "Bearer " + token);
            }

            // Post  content
            outputStream = connection.getOutputStream();
            outputStream.write(bytes);
            outputStream.close();

            // Read the response
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(connection.getInputStream())));
            StringBuffer sb = new StringBuffer();
            String str = br.readLine();
            while (str != null) {
                sb.append(str);
                str = br.readLine();
            }
            br.close();

            return sb.toString();
        } catch (IOException e) {

            int rs = connection.getResponseCode();
            if (rs >= 400 && rs < 500) {
                displayAuthorizationDialog(url.getHost());
            } else {
                MessageUtils.showErrorMessage("Error accessing resource", e);
                e.printStackTrace();
            }

            return null;
        }
    }

    static void displayAuthorizationDialog(String host) {

        String message = "The requested resource at '"  + host + "' requires authorization.";
        Icon icon = null;
        int option = JOptionPane.showOptionDialog(IGV.getMainFrame(),
                message,
                "Error",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                icon,
                new String[]{"Cancel", "Authorize"},
                JOptionPane.YES_OPTION
        );

        if(option == 1) {
            try {
                OAuthUtils.getInstance().fetchAuthCode();
            } catch (Exception e) {
                MessageUtils.showErrorMessage("Error fetching oAuth token", e);
                log.error("Error fetching oAuth tokens", e);
            }
        }

    }


    static Map<String, String> genomeIdMap = new HashMap<String, String>();  // A hack until readset meta data is available

    static {
        genomeIdMap = new HashMap<String, String>();
        genomeIdMap.put("Google 10473108253681171589", "hg19");
        genomeIdMap.put("Google 383928317087", "hg19");
        genomeIdMap.put("Google 461916304629", "hg19");
        genomeIdMap.put("Google 337315832689", "hg19");

        genomeIdMap.put("NCBI SRP034507", "M74568");
        genomeIdMap.put("NCBI SRP029392", "NC_004917");

    }


}
