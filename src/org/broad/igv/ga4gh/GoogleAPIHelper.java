package org.broad.igv.ga4gh;

import com.google.gson.*;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.Pair;

import java.awt.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.*;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Helper class for  Google API prototype
 * <p/>
 * Created by jrobinso on 8/15/14.
 */
public class GoogleAPIHelper {
    public static final String RESOURCE_TYPE = "ga4gh";

    // Magic dataset id (1000 genomes)

    final static String datasetId = "376902546192";

    final static Map<String, List<Pair<String, String>>> readsetCache = new HashMap<String, List<Pair<String, String>>>();

//
//    public static void main(String[] args) throws IOException {
//        List<Pair<String, String>> idNamePairs = readsetSearch(datasetId);
//        (new GoogleAPILoadDialog(null, idNamePairs)).setVisible(true);
//    }


    public static void openLoadDialog(IGV igv, Frame frame) throws IOException {

        List<Pair<String, String>> idNamePairs = readsetSearch(datasetId);
        GoogleAPILoadDialog dlg = (new GoogleAPILoadDialog(frame, idNamePairs));
        dlg.setModal(true);
        dlg.setVisible(true);
        dlg.dispose();
    }

    public static List<Pair<String, String>> readsetSearch(String datasetId) throws IOException {

        List<Pair<String, String>> idNamePairs = readsetCache.get(datasetId);
        if (idNamePairs == null) {

            idNamePairs = new ArrayList();

            // Loop through pages
            int maxPages = 100;
            JsonPrimitive pageToken = null;
            while (maxPages-- > 0) {
                String contentToPost = "{" +
                        "\"datasetIds\": [" + datasetId + "]" +
                        (pageToken == null ? "" : ", \"pageToken\": " + pageToken) +
                        ", \"maxResults\": 256" +
                        "}";

                String result = doPost("/readsets/search", contentToPost, "&fields=readsets(id,name),nextPageToken");

                JsonParser parser = new JsonParser();
                JsonObject obj = parser.parse(result).getAsJsonObject();
                JsonArray readsets = obj.getAsJsonArray("readsets");

                Iterator<JsonElement> iter = readsets.iterator();

                while (iter.hasNext()) {
                    JsonElement next = iter.next();
                    JsonObject jobj = next.getAsJsonObject();
                    idNamePairs.add(new Pair(jobj.get("id").getAsString(), jobj.get("name").getAsString()));
                }

                pageToken = obj.getAsJsonPrimitive("nextPageToken");
                if (pageToken == null) break;
            }

            Collections.sort(idNamePairs, new Comparator<Pair<String, String>>() {
                @Override
                public int compare(Pair<String, String> o1, Pair<String, String> o2) {
                    return o1.getSecond().compareTo(o2.getSecond());
                }
            });

            readsetCache.put(datasetId, idNamePairs);
        }

        return idNamePairs;
    }


    public static List<Alignment> reads(String readsetId, String chr, int start, int end) throws IOException {

        List<Alignment> alignments = new ArrayList<Alignment>(10000);
        int maxPages = 10000;
        JsonPrimitive pageToken = null;
        StringBuffer result = new StringBuffer();

        while (maxPages-- > 0) {
            String contentToPost = "{" +
                    "readsetIds: [\"" + readsetId + "\"]" +
                    ", sequenceName: " + chr +
                    ", sequenceStart: " + start +
                    ", sequenceEnd: " + end +
                    (pageToken == null ? "" : ", \"pageToken\": " + pageToken) +
                    "}";

            String readString =  doPost("/reads/search", contentToPost, "");

            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(readString).getAsJsonObject();

            JsonArray reads = obj.getAsJsonArray("reads");

            Iterator<JsonElement> iter = reads.iterator();
            while (iter.hasNext()) {
                JsonElement next = iter.next();
                Ga4ghAlignment alignment = new Ga4ghAlignment(next.getAsJsonObject());
                alignments.add(alignment);
            }

            pageToken = obj.getAsJsonPrimitive("nextPageToken");
            if (pageToken == null) break;

        }

        System.out.println("# pages= " + (10000 - maxPages));

        return alignments;

    }

    private static String doPost(String command, String content, String fields) throws IOException {

        String authKey = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_API_KEY);
        String baseURL = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_BASE_URL);
        URL url = new URL(baseURL + command + "?key=" + authKey + fields);


        byte[] bytes = content.getBytes();

        // Create a URLConnection
        HttpURLConnection connection = (HttpURLConnection) url.openConnection();
        connection.setUseCaches(false);
        connection.setDoInput(true);
        connection.setDoOutput(true);
        connection.setRequestMethod("POST");
        //connection.setRequestProperty("Content-Length", "" + bytes.length);
        connection.setRequestProperty("Content-Type", "application/json");
        connection.setRequestProperty("Cache-Control", "no-cache");
        connection.setRequestProperty("Accept-Encoding", "gzip");
        connection.setRequestProperty("User-Agent", "IGV (gzip)");

        // Post  content
        java.io.OutputStream stream = connection.getOutputStream();
        stream.write(bytes);
        stream.close();

        // Read the response
        java.io.BufferedReader br = new java.io.BufferedReader(new InputStreamReader(new GZIPInputStream(connection.getInputStream())));
        StringBuffer sb = new StringBuffer();
        String str = br.readLine();
        while (str != null) {
            sb.append(str);
            str = br.readLine();
        }
        br.close();

        return sb.toString();
    }


}
