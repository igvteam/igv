package org.igv.session;

import org.apache.commons.io.IOUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.util.ResourceLocator;
import org.json.JSONObject;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

public class JSONSessionReader implements SessionReader{

    private final IGV igv;

    public JSONSessionReader(IGV igv) {
        this.igv = igv;
    }

    @Override
    public void loadSession(InputStream inputStream, Session session, String sessionPath) throws IOException {

        String jsonString = IOUtils.toString(inputStream, StandardCharsets.UTF_8);
        JSONObject jsonObject = new JSONObject(jsonString);

        if(jsonObject.has("genome")) {
            GenomeManager.getInstance().loadGenomeById(jsonObject.getString("genome"));
        }

        if(jsonObject.has("locus")) {
            igv.goToLocus(jsonObject.getString("locus"));
        }

        if (jsonObject.has("tracks")) {

            List<CompletableFuture<List<Track>>> futures = new ArrayList<>();
            List<JSONObject> trackJsons = new ArrayList<>();

            jsonObject.getJSONArray("tracks").forEach(trackObj -> {
                JSONObject trackJson = (JSONObject) trackObj;
                if (trackJson.has("url")) {
                    trackJsons.add(trackJson);
                    String url = trackJson.getString("url");
                    ResourceLocator locator = new ResourceLocator(url);
                    futures.add(CompletableFuture.supplyAsync(() -> igv.load(locator)));
                }
            });

            CompletableFuture.allOf(futures.toArray(new CompletableFuture[0])).join();

            for (int i = 0; i < futures.size(); i++) {
                try {
                    List<Track> tracks = futures.get(i).get();
                    JSONObject trackJson = trackJsons.get(i);
                    for (Track t : tracks) {
                        t.unmarshalJSON(trackJson);
                        igv.addTrack(t);
                    }
                } catch (Exception e) {
                    // Log error or handle exception
                }
            }
        }

    }
}
