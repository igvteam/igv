package org.igv.session;

import org.apache.commons.io.IOUtils;
import org.igv.data.CombinedDataSource;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.CombinedDataTrack;
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

                } else if (trackJson.has("type") && "combined".equals(trackJson.getString("type"))) {
                    // Combined track - defer creation until referenced tracks are loaded
                    topLevelDescriptors.add(new TrackDescriptor(trackJson));

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

            // Map to store id -> track for combined track creation
            Map<String, Track> trackById = new HashMap<>();

            // Track the current index for insertion and remember combined track positions
            int currentIndex = 0;
            List<Integer> combinedTrackIndices = new ArrayList<>();
            List<TrackDescriptor> combinedDescriptors = new ArrayList<>();

            // First pass: add non-combined tracks in order, record combined track positions
            for (TrackDescriptor descriptor : topLevelDescriptors) {
                try {
                    if (descriptor.isCombined()) {
                        // Remember the position for this combined track
                        combinedTrackIndices.add(currentIndex);
                        combinedDescriptors.add(descriptor);
                        currentIndex++;

                    } else if (descriptor.isMerged()) {
                        // Assemble merged track from loaded child tracks
                        Collection<DataTrack> childTracks = new ArrayList<>();
                        for (int futureIndex : descriptor.childFutureIndices) {
                            List<Track> loadedTracks = allTrackFutures.get(futureIndex).get();
                            for (Track track : loadedTracks) {
                                if (track instanceof DataTrack) {
                                    childTracks.add((DataTrack) track);
                                    trackById.put(track.getId(), track);
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
                            trackById.put(mergedTracks.getId(), mergedTracks);
                            currentIndex++;
                        } else {
                            log.warn("Merged track has no valid child tracks: " + descriptor.trackJson.optString("name", "unknown"));
                        }

                    } else {
                        // Add non-merged track(s)
                        List<Track> tracks = allTrackFutures.get(descriptor.singleFutureIndex).get();
                        for (Track t : tracks) {
                            igv.addTrack(t);
                            trackById.put(t.getId(), t);
                            currentIndex++;
                        }
                    }
                } catch (Exception e) {
                    log.error("Error processing track descriptor", e);
                }
            }

            // Second pass: create and insert combined tracks at their original positions
            for (int i = 0; i < combinedDescriptors.size(); i++) {
                TrackDescriptor descriptor = combinedDescriptors.get(i);
                int insertIndex = combinedTrackIndices.get(i);
                try {
                    CombinedDataTrack combinedTrack = createCombinedDataTrack(descriptor.trackJson, trackById);
                    if (combinedTrack != null) {
                        combinedTrack.unmarshalJSON(descriptor.trackJson);
                        igv.addTracks(List.of(combinedTrack), insertIndex);
                    }
                } catch (Exception e) {
                    log.error("Error creating combined track", e);
                }
            }
        }
    }

    /**
     * Create a CombinedDataTrack from the JSON definition and the map of loaded tracks.
     *
     * @param trackJson   the JSON object containing track1, track2, and op properties
     * @param trackById   map of track id to Track object
     * @return the created CombinedDataTrack, or null if referenced tracks are not found
     */
    private CombinedDataTrack createCombinedDataTrack(JSONObject trackJson, Map<String, Track> trackById) {
        String track1Id = trackJson.getString("track1");
        String track2Id = trackJson.getString("track2");
        String opString = trackJson.getString("op");

        Track track1 = trackById.get(track1Id);
        Track track2 = trackById.get(track2Id);

        if (track1 == null) {
            log.error("Combined track references unknown track1: " + track1Id);
            return null;
        }
        if (track2 == null) {
            log.error("Combined track references unknown track2: " + track2Id);
            return null;
        }
        if (!(track1 instanceof DataTrack)) {
            log.error("Combined track track1 is not a DataTrack: " + track1.getClass().getName());
            return null;
        }
        if (!(track2 instanceof DataTrack)) {
            log.error("Combined track track2 is not a DataTrack: " + track2.getClass().getName());
            return null;
        }

        CombinedDataSource.Operation operation = CombinedDataSource.Operation.valueOf(opString);
        CombinedDataSource dataSource = new CombinedDataSource((DataTrack) track1, (DataTrack) track2, operation);

        String id = trackJson.optString("id", UUID.randomUUID().toString());
        String name = trackJson.optString("name", "Combined Track");

        return new CombinedDataTrack(dataSource, id, name);
    }

    /**
     * Descriptor for a top-level track (merged, combined, or non-merged) that maintains
     * references to the futures that will provide the loaded tracks.
     */
    private static class TrackDescriptor {
        enum Type { SINGLE, MERGED, COMBINED }

        final Type type;
        final JSONObject trackJson;
        final List<Integer> childFutureIndices; // For merged tracks
        final List<JSONObject> childJsons;      // For merged tracks
        final int singleFutureIndex;            // For non-merged tracks

        // Constructor for merged track
        TrackDescriptor(JSONObject trackJson, List<Integer> childFutureIndices, List<JSONObject> childJsons) {
            this.type = Type.MERGED;
            this.trackJson = trackJson;
            this.childFutureIndices = childFutureIndices;
            this.childJsons = childJsons;
            this.singleFutureIndex = -1;
        }

        // Constructor for non-merged track
        TrackDescriptor(JSONObject trackJson, int futureIndex) {
            this.type = Type.SINGLE;
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = futureIndex;
        }

        // Constructor for combined track (no futures, created after other tracks)
        TrackDescriptor(JSONObject trackJson) {
            this.type = Type.COMBINED;
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = -1;
        }

        boolean isMerged() {
            return type == Type.MERGED;
        }

        boolean isCombined() {
            return type == Type.COMBINED;
        }
    }
}
