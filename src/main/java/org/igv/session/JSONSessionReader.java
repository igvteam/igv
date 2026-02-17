package org.igv.session;

import org.apache.commons.io.IOUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.DataTrack;
import org.igv.track.MergedTracks;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.util.ResourceLocator;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.CompletableFuture;

public class JSONSessionReader implements SessionReader {

    private static Logger log = LogManager.getLogger(JSONSessionReader.class);

    private final IGV igv;

    public JSONSessionReader(IGV igv) {
        this.igv = igv;
    }

    @Override
    public void loadSession(InputStream inputStream, Session session, String sessionPath) throws IOException {

        String jsonString = IOUtils.toString(inputStream, StandardCharsets.UTF_8);
        JSONObject jsonObject = new JSONObject(jsonString);

        if (jsonObject.has("genome")) {
            GenomeManager.getInstance().loadGenomeById(jsonObject.getString("genome"));
        }

        if (jsonObject.has("locus")) {
            igv.goToLocus(jsonObject.getString("locus"));
        }

        if (jsonObject.has("tracks")) {

            JSONArray tracksArray = jsonObject.getJSONArray("tracks");

            // Collect all individual track load futures (both top-level and merged children)
            List<CompletableFuture<List<Track>>> allTrackFutures = new ArrayList<>();

            // Track descriptors maintain the order and structure of top-level tracks
            List<TrackDescriptor> topLevelDescriptors = new ArrayList<>();

            for (int i = 0; i < tracksArray.length(); i++) {
                JSONObject trackJson = tracksArray.getJSONObject(i);

                if (trackJson.has("type") && "merged".equals(trackJson.getString("type"))) {
                    // Merged track - collect child track futures
                    List<Integer> childFutureIndices = new ArrayList<>();
                    List<JSONObject> childJsons = new ArrayList<>();

                    JSONArray childTracksArray = trackJson.getJSONArray("tracks");
                    for (int j = 0; j < childTracksArray.length(); j++) {
                        JSONObject childTrackJson = childTracksArray.getJSONObject(j);
                        if (childTrackJson.has("url")) {
                            String url = childTrackJson.getString("url");
                            ResourceLocator locator = new ResourceLocator(url);

                            int futureIndex = allTrackFutures.size();
                            childFutureIndices.add(futureIndex);
                            childJsons.add(childTrackJson);

                            allTrackFutures.add(CompletableFuture.supplyAsync(() -> {
                                List<Track> tracks = igv.load(locator);
                                tracks.forEach(track -> track.unmarshalJSON(childTrackJson));
                                return tracks;
                            }).exceptionally(e -> {
                                log.error("Error loading child track: " + url, e);
                                return Collections.emptyList();
                            }));
                        }
                    }

                    topLevelDescriptors.add(new TrackDescriptor(trackJson, childFutureIndices, childJsons));

                } else if (trackJson.has("url")) {
                    // Non-merged top-level track
                    String url = trackJson.getString("url");
                    ResourceLocator locator = new ResourceLocator(url);

                    int futureIndex = allTrackFutures.size();
                    topLevelDescriptors.add(new TrackDescriptor(trackJson, futureIndex));

                    allTrackFutures.add(CompletableFuture.supplyAsync(() -> {
                        List<Track> tracks = igv.load(locator);
                        for (Track t : tracks) {
                            t.unmarshalJSON(trackJson);
                        }
                        return tracks;
                    }).exceptionally(e -> {
                        log.error("Error loading track: " + url, e);
                        return Collections.emptyList();
                    }));
                }
            }

            // Wait for all track loads to complete in parallel (use exceptionally handlers above to prevent failures)
            try {
                CompletableFuture.allOf(allTrackFutures.toArray(new CompletableFuture[0])).join();
            } catch (Exception e) {
                log.error("Error waiting for track futures to complete", e);
            }

            // Reconstruct top-level tracks in order
            for (TrackDescriptor descriptor : topLevelDescriptors) {
                try {
                    if (descriptor.isMerged()) {
                        // Assemble merged track from loaded child tracks
                        Collection<DataTrack> childTracks = new ArrayList<>();
                        for (int futureIndex : descriptor.childFutureIndices) {
                            List<Track> loadedTracks = allTrackFutures.get(futureIndex).get();
                            for (Track track : loadedTracks) {
                                if (track instanceof DataTrack) {
                                    childTracks.add((DataTrack) track);
                                } else {
                                    log.warn("Expected DataTrack but got: " + track.getClass().getName());
                                }
                            }
                        }

                        if (!childTracks.isEmpty()) {
                            String id = descriptor.trackJson.optString("id", UUID.randomUUID().toString());
                            String name = descriptor.trackJson.optString("name", "Merged Track");
                            MergedTracks mergedTracks = new MergedTracks(id, name, childTracks);
                            mergedTracks.unmarshalJSON(descriptor.trackJson);
                            igv.addTrack(mergedTracks);
                        } else {
                            log.warn("Merged track has no valid child tracks: " + descriptor.trackJson.optString("name", "unknown"));
                        }

                    } else {
                        // Add non-merged track(s)
                        List<Track> tracks = allTrackFutures.get(descriptor.singleFutureIndex).get();
                        for (Track t : tracks) {
                            igv.addTrack(t);
                        }
                    }
                } catch (Exception e) {
                    log.error("Error processing track descriptor", e);
                }
            }
        }
    }

    /**
     * Descriptor for a top-level track (merged or non-merged) that maintains
     * references to the futures that will provide the loaded tracks.
     */
    private static class TrackDescriptor {
        final JSONObject trackJson;
        final List<Integer> childFutureIndices; // For merged tracks
        final List<JSONObject> childJsons;      // For merged tracks
        final int singleFutureIndex;            // For non-merged tracks

        // Constructor for merged track
        TrackDescriptor(JSONObject trackJson, List<Integer> childFutureIndices, List<JSONObject> childJsons) {
            this.trackJson = trackJson;
            this.childFutureIndices = childFutureIndices;
            this.childJsons = childJsons;
            this.singleFutureIndex = -1;
        }

        // Constructor for non-merged track
        TrackDescriptor(JSONObject trackJson, int futureIndex) {
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = futureIndex;
        }

        boolean isMerged() {
            return childFutureIndices != null;
        }
    }
}
