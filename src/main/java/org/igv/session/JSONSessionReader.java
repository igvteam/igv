package org.igv.session;

import org.apache.commons.io.IOUtils;
import org.igv.data.CombinedDataSource;
import org.igv.feature.RegionOfInterest;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.*;
import org.igv.ui.IGV;
import org.igv.util.FilterElement;
import org.igv.util.ResourceLocator;
import org.igv.util.TrackFilter;
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
            if (!GenomeManager.getInstance().getCurrentGenome().getId().equals(jsonObject.getString("genome"))) {
                GenomeManager.getInstance().loadGenomeById(jsonObject.getString("genome"));
            }

        } else if (jsonObject.has("reference")) {
            JSONObject genomeJson = jsonObject.getJSONObject("reference");
            GenomeConfig genomeConfig = GenomeConfig.getGenomeConfig(genomeJson);
            Genome newGenome = new Genome(genomeConfig);
            if (IGV.hasInstance()) {
                IGV.getInstance().resetSession(null);
            }
            GenomeManager.getInstance().setCurrentGenome(newGenome);

        } else {
            throw new RuntimeException("Session file is missing required 'genome' or 'reference' property");
        }

        // Cache genome tracks if provided.  Genome tracks will always include sequence, but may also include
        // annotation tracks.  We remove them here and re-add after loading all other tracks to maintain the specified
        // track order.
        Map<String, Track> genomeTrackMap = new HashMap<>();
        List<Track> genomeTracks = igv.getAllTracks();
        for (Track track : genomeTracks) {
            genomeTrackMap.put(getId(track), track);
        }
        igv.clearTrackPanels();

        if (jsonObject.has("locus")) {
            Object locusValue = jsonObject.get("locus");
            if (locusValue instanceof JSONArray) {
                // Array of locus values
                JSONArray locusArray = (JSONArray) locusValue;
                List<String> loci = new ArrayList<>();
                for (int i = 0; i < locusArray.length(); i++) {
                    loci.add(locusArray.getString(i));
                }
                igv.goToLocus(String.join(" ", loci));
            } else {
                // Single locus string (may be space-delimited string of multiple loci)
                igv.goToLocus(locusValue.toString());
            }
        }

        if(jsonObject.has("groupSampleAttribute")) {
            session.setGroupByAttribute(jsonObject.getString("groupSampleAttribute"));
        } else if(jsonObject.has("groupTracksBy")) {
            session.setGroupByAttribute(jsonObject.getString("groupTracksBy"));
        }

        if (jsonObject.has("tracks")) {

            JSONArray tracksArray = jsonObject.getJSONArray("tracks");

            // Collect all individual track load futures (both top-level and merged children)
            List<CompletableFuture<List<Track>>> allTrackFutures = new ArrayList<>();

            // Map URL to future index to avoid loading the same URL multiple times
            Map<String, Integer> urlToFutureIndex = new HashMap<>();

            // Track descriptors maintain the order and structure of top-level tracks
            List<TrackDescriptor> topLevelDescriptors = new ArrayList<>();

            for (int i = 0; i < tracksArray.length(); i++) {

                JSONObject trackJson = tracksArray.getJSONObject(i);
                String type = trackJson.optString("type", null);
                String id = getId(trackJson);
                if (id != null && genomeTrackMap.containsKey(id)) {

                    // This track is already loaded as a genome track - create descriptor for correct positioning
                    Track existingTrack = genomeTrackMap.remove(id);
                    topLevelDescriptors.add(new TrackDescriptor(trackJson, existingTrack, i));
                    continue;
                }

                if ("merged".equals(type)) {
                    // Merged track - collect child track futures
                    List<Integer> childFutureIndices = new ArrayList<>();
                    List<JSONObject> childJsons = new ArrayList<>();

                    JSONArray childTracksArray = trackJson.getJSONArray("tracks");
                    for (int j = 0; j < childTracksArray.length(); j++) {
                        JSONObject childTrackJson = childTracksArray.getJSONObject(j);
                        if (childTrackJson.has("url")) {
                            String url = childTrackJson.getString("url");

                            int futureIndex;
                            if (urlToFutureIndex.containsKey(url)) {
                                futureIndex = urlToFutureIndex.get(url);
                            } else {
                                ResourceLocator locator = new ResourceLocator(url);
                                String format = childTrackJson.optString("format", null);
                                if (format != null) {
                                    locator.setFormat(format);
                                }
                                futureIndex = allTrackFutures.size();
                                urlToFutureIndex.put(url, futureIndex);
                                allTrackFutures.add(CompletableFuture.supplyAsync(() -> igv.load(locator))
                                        .exceptionally(e -> {
                                            log.error("Error loading child track: " + url, e);
                                            return Collections.emptyList();
                                        }));
                            }

                            childFutureIndices.add(futureIndex);
                            childJsons.add(childTrackJson);
                        }
                    }

                    topLevelDescriptors.add(new TrackDescriptor(trackJson, childFutureIndices, childJsons, i));

                } else if ("combined".equals(type)) {
                    // Combined track - defer creation until referenced tracks are loaded
                    topLevelDescriptors.add(new TrackDescriptor(trackJson, i, TrackDescriptor.Type.COMBINED));

                } else if("blat".equals(type)) {
                    topLevelDescriptors.add(TrackDescriptor.createBlatDescriptor(trackJson, i));

                } else if("sequence".equals(type)) {
                    // Sequence track - find in genomeTrackMap by type
                    Track sequenceTrack = null;
                    for (Map.Entry<String, Track> entry : genomeTrackMap.entrySet()) {
                        if (entry.getValue() instanceof SequenceTrack) {
                            sequenceTrack = entry.getValue();
                            genomeTrackMap.remove(entry.getKey());
                            break;
                        }
                    }
                    if (sequenceTrack != null) {
                        topLevelDescriptors.add(new TrackDescriptor(trackJson, sequenceTrack, i));
                    } else {
                        log.warn("Sequence track not found in genome tracks");
                    }
                } else if("motif".equals(type)) {
                    topLevelDescriptors.add(TrackDescriptor.createMotifDescriptor(trackJson, i));

                } else if (trackJson.has("url")) {
                    // Non-merged top-level track
                    String url = trackJson.getString("url");

                    int futureIndex;
                    if (urlToFutureIndex.containsKey(url)) {
                        // Reuse existing future for this URL
                        futureIndex = urlToFutureIndex.get(url);
                    } else {
                        // Create new future for this URL
                        ResourceLocator locator = new ResourceLocator(url);
                        String format = trackJson.optString("format", null);
                        if (format != null) {
                            locator.setFormat(format);
                        }
                        futureIndex = allTrackFutures.size();
                        urlToFutureIndex.put(url, futureIndex);
                        allTrackFutures.add(CompletableFuture.supplyAsync(() -> igv.load(locator))
                                .exceptionally(e -> {
                                    log.error("Error loading track: " + url, e);
                                    return Collections.emptyList();
                                }));
                    }

                    topLevelDescriptors.add(new TrackDescriptor(trackJson, futureIndex, i));
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

            // Track combined descriptors for second pass (combined tracks depend on other tracks being loaded first)
            List<TrackDescriptor> combinedDescriptors = new ArrayList<>();

            // First pass: add all non-combined tracks, setting order property from descriptor
            for (TrackDescriptor descriptor : topLevelDescriptors) {
                try {
                    if (descriptor.isCombined()) {
                        // Defer combined tracks to second pass
                        combinedDescriptors.add(descriptor);

                    } else if (descriptor.isBlat()) {
                        // Blat track - create, set order, and add
                        BlatTrack blatTrack = new BlatTrack();
                        blatTrack.setOrder(descriptor.order);
                        blatTrack.unmarshalJSON(descriptor.trackJson);
                        igv.addTrack(blatTrack);
                        trackById.put(getId(blatTrack), blatTrack);

                    } else if (descriptor.isMotif()) {
                        // Motif track - create, set order, and add
                        MotifTrack motifTrack = new MotifTrack();
                        motifTrack.setOrder(descriptor.order);
                        motifTrack.unmarshalJSON(descriptor.trackJson);
                        igv.addTrack(motifTrack);
                        trackById.put(getId(motifTrack), motifTrack);

                    } else if (descriptor.isGenome()) {
                        // Genome track - set order, unmarshal, and add
                        Track t = descriptor.genomeTrack;
                        t.setOrder(descriptor.order);
                        t.unmarshalJSON(descriptor.trackJson);
                        igv.addTrack(t);
                        trackById.put(getId(t), t);

                    } else if (descriptor.isMerged()) {
                        // Assemble merged track from loaded child tracks
                        Collection<DataTrack> childTracks = new ArrayList<>();
                        for (int j = 0; j < descriptor.childFutureIndices.size(); j++) {
                            int futureIndex = descriptor.childFutureIndices.get(j);
                            JSONObject childJson = descriptor.childJsons.get(j);

                            List<Track> loadedTracks = allTrackFutures.get(futureIndex).get();
                            for (Track track : loadedTracks) {
                                if (track instanceof DataTrack) {
                                    track.unmarshalJSON(childJson);
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
                            mergedTracks.setOrder(descriptor.order);
                            mergedTracks.unmarshalJSON(descriptor.trackJson);
                            igv.addTrack(mergedTracks);
                            trackById.put(mergedTracks.getId(), mergedTracks);
                        } else {
                            log.warn("Merged track has no valid child tracks: " + descriptor.trackJson.optString("name", "unknown"));
                        }

                    } else {
                        // Non-merged track(s) - match loaded tracks to trackJson using name, id, or type
                        int futureIndex = descriptor.singleFutureIndex;
                        List<Track> tracks = allTrackFutures.get(futureIndex).get();

                        // Find the track that matches this descriptor's trackJson
                        String jsonName = descriptor.trackJson.optString("name", null);
                        String jsonId = descriptor.trackJson.optString("id", null);
                        String jsonType = descriptor.trackJson.optString("type", null);
                        String jsonFormat = descriptor.trackJson.optString("format", null);

                        for (Track t : tracks) {
                            // Match track to trackJson using name, id, type, or format
                            boolean matches = (jsonName != null && jsonName.equals(t.getName())) ||
                                    (jsonId != null && jsonId.equals(t.getId())) ||
                                    (jsonType != null && jsonType.equalsIgnoreCase(t.getType().toString())) ||
                                    (jsonFormat != null && jsonFormat.equalsIgnoreCase(t.getType().toString())) ||
                                    (tracks.size() == 1);

                            // Use a unique key that includes the track type to allow multiple tracks from same URL
                            String trackKey = getId(t) + "|" + t.getType();

                            if (matches && !trackById.containsKey(trackKey)) {
                                t.unmarshalJSON(descriptor.trackJson);
                                t.setOrder(descriptor.order);
                                igv.addTrack(t);
                                trackById.put(trackKey, t);
                                break; // Only match one track per descriptor
                            }
                        }
                    }
                } catch (Exception e) {
                    log.error("Error processing track descriptor", e);
                }
            }

            // Second pass: create and add combined tracks
            for (TrackDescriptor descriptor : combinedDescriptors) {
                try {
                    CombinedDataTrack combinedTrack = createCombinedDataTrack(descriptor.trackJson, trackById);
                    if (combinedTrack != null) {
                        combinedTrack.setOrder(descriptor.order);
                        combinedTrack.unmarshalJSON(descriptor.trackJson);
                        igv.addTrack(combinedTrack);
                    }
                } catch (Exception e) {
                    log.error("Error creating combined track", e);
                }
            }

            // Add any remaining tracks that are explicitly defined in the genome or reference section, but were not
            // referenced in the tracks section.
            for (Track track : genomeTrackMap.values()) {
                igv.addTrack(track);
            }
        }

        if(jsonObject.has("sampleinfo")) {
            JSONArray sampleInfoArray = jsonObject.getJSONArray("sampleinfo");
            for (int i = 0; i < sampleInfoArray.length(); i++) {
                JSONObject sampleInfoJson = sampleInfoArray.getJSONObject(i);
                if (sampleInfoJson.has("url")) {
                    String url = sampleInfoJson.getString("url");
                    ResourceLocator locator = new ResourceLocator(url);
                    try {
                        AttributeManager.getInstance().loadSampleInfo(locator);
                    } catch (Exception e) {
                        log.error("Error loading sample info: " + url, e);
                    }
                }
            }
        }

        if (jsonObject.has("roi")) {
            JSONArray roiArray = jsonObject.getJSONArray("roi");
            for (int i = 0; i < roiArray.length(); i++) {
                JSONObject roiObject = roiArray.getJSONObject(i);
                if (roiObject.has("features")) {
                    JSONArray featuresArray = roiObject.getJSONArray("features");
                    for (int j = 0; j < featuresArray.length(); j++) {
                        JSONObject featureJson = featuresArray.getJSONObject(j);
                        String chr = featureJson.getString("chr");
                        int start = featureJson.getInt("start");
                        int end = featureJson.getInt("end");
                        String description = featureJson.optString("description", null);
                        if(description == null) {
                            // "name" is used for description in igv.js sessions
                            description = featureJson.optString("name", null);
                        }
                        RegionOfInterest roi = new RegionOfInterest(chr, start, end, description);
                        session.addRegionOfInterest(roi);
                    }
                }
            }
        }

        // Restore sort attributes after all tracks and sample info are loaded
        if (jsonObject.has("sortSamples")) {
            JSONObject sortJson = jsonObject.getJSONObject("sortSamples");
            if (sortJson.has("attributeNames") && sortJson.has("ascending")) {
                JSONArray attributeNamesArray = sortJson.getJSONArray("attributeNames");
                JSONArray ascendingArray = sortJson.getJSONArray("ascending");

                String[] attributeNames = new String[attributeNamesArray.length()];
                boolean[] ascending = new boolean[ascendingArray.length()];

                for (int i = 0; i < attributeNamesArray.length(); i++) {
                    attributeNames[i] = attributeNamesArray.getString(i);
                    ascending[i] = ascendingArray.getBoolean(i);
                }

                // Apply the sorting
                IGV.getInstance().sortAllTracksByAttributes(attributeNames, ascending);
            }
        }

        // Restore filter attributes after all tracks and sample info are loaded
        if (jsonObject.has("filterSamples")) {
            JSONObject filterJson = jsonObject.getJSONObject("filterSamples");

            // Determine match mode (all or any)
            String matchMode = filterJson.optString(SessionAttribute.FILTER_MATCH, "all");
            boolean matchAll = "all".equalsIgnoreCase(matchMode);

            // Read filter elements
            List<FilterElement> filterElements = new ArrayList<>();
            if (filterJson.has("elements")) {
                JSONArray elementsArray = filterJson.getJSONArray("elements");
                for (int i = 0; i < elementsArray.length(); i++) {
                    JSONObject elementJson = elementsArray.getJSONObject(i);

                    String attributeKey = elementJson.getString(SessionAttribute.ITEM);
                    String operatorString = elementJson.getString(SessionAttribute.OPERATOR);
                    String value = elementJson.getString(SessionAttribute.VALUE);

                    // Convert operator string to Operator enum
                    FilterElement.Operator operator = getOperatorFromValue(operatorString);
                    if (operator != null) {
                        filterElements.add(new FilterElement(attributeKey, operator, value));
                    } else {
                        log.warn("Unknown filter operator: " + operatorString);
                    }
                }
            }

            if (!filterElements.isEmpty()) {
                TrackFilter trackFilter = new TrackFilter(matchAll, filterElements);
                session.setFilter(trackFilter);

                // Apply the filter to all tracks
                for (Track track : IGV.getInstance().getAllTracks()) {
              //      track.groupSamplesByAttribute();
                }
            }
        }
    }

    /**
     * Convert an operator value string to the corresponding Operator enum.
     */
    private static FilterElement.Operator getOperatorFromValue(String value) {
        for (FilterElement.Operator op : FilterElement.Operator.values()) {
            if (op.getValue().equalsIgnoreCase(value)) {
                return op;
            }
        }
        return null;
    }


    /**
     * Return an identifier for a track JSON object. The "id" property is not reliable, so we check "url" and "path"
     * properties, which are commonly used.
     *
     * @param trackJson
     * @return
     */
    private static String getId(JSONObject trackJson) {
        return trackJson.has("url") ? trackJson.getString("url") :
                trackJson.has("path") ? trackJson.getString("path") :
                trackJson.optString("id", null);
    }

    /**
     * Return an identifier for a loaded Track.  We check the resource locator path first, as this is more likely
     * to be consistent with the JSON definition than the track's internal id (which may be a generated UUID).
     *
     * @param track
     * @return
     */
    private static String getId(Track track) {
        ResourceLocator locator = track.getResourceLocator();
        return locator != null && locator.getPath() != null ? locator.getPath() : track.getId();
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
     * Descriptor for a top-level track (merged, combined, blat, genome, or non-merged) that maintains
     * references to the futures that will provide the loaded tracks.
     */
    private static class TrackDescriptor {
        enum Type {SINGLE, MERGED, COMBINED, BLAT, MOTIF, GENOME}

        final Type type;
        final JSONObject trackJson;
        final List<Integer> childFutureIndices; // For merged tracks
        final List<JSONObject> childJsons;      // For merged tracks
        final int singleFutureIndex;            // For non-merged tracks
        final Track genomeTrack;                // For genome tracks
        final int fileIndex;                    // Original index in the JSON file
        final long order;                       // Explicit order value (0 if not specified)

        // Constructor for merged track
        TrackDescriptor(JSONObject trackJson, List<Integer> childFutureIndices, List<JSONObject> childJsons, int fileIndex) {
            this.type = Type.MERGED;
            this.trackJson = trackJson;
            this.childFutureIndices = childFutureIndices;
            this.childJsons = childJsons;
            this.singleFutureIndex = -1;
            this.genomeTrack = null;
            this.fileIndex = fileIndex;
            this.order = trackJson.optLong("order", 0);
        }

        // Constructor for non-merged track
        TrackDescriptor(JSONObject trackJson, int futureIndex, int fileIndex) {
            this.type = Type.SINGLE;
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = futureIndex;
            this.genomeTrack = null;
            this.fileIndex = fileIndex;
            this.order = trackJson.optLong("order", 0);
        }

        // Constructor for combined track (no futures, created after other tracks)
        TrackDescriptor(JSONObject trackJson, int fileIndex, Type type) {
            this.type = type;
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = -1;
            this.genomeTrack = null;
            this.fileIndex = fileIndex;
            this.order = trackJson.optLong("order", 0);
        }

        // Constructor for blat track (created inline, no url loading needed)
        static TrackDescriptor createBlatDescriptor(JSONObject trackJson, int fileIndex) {
            return new TrackDescriptor(Type.BLAT, trackJson, fileIndex);
        }

        // Constructor for motif track (created inline, no url loading needed)
        static TrackDescriptor createMotifDescriptor(JSONObject trackJson, int fileIndex) {
            return new TrackDescriptor(Type.MOTIF, trackJson, fileIndex);
        }

        // Private constructor for special track types
        private TrackDescriptor(Type type, JSONObject trackJson, int fileIndex) {
            this.type = type;
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = -1;
            this.genomeTrack = null;
            this.fileIndex = fileIndex;
            this.order = trackJson.optLong("order", 0);
        }

        // Constructor for genome track (already loaded, just needs positioning)
        TrackDescriptor(JSONObject trackJson, Track genomeTrack, int fileIndex) {
            this.type = Type.GENOME;
            this.trackJson = trackJson;
            this.childFutureIndices = null;
            this.childJsons = null;
            this.singleFutureIndex = -1;
            this.genomeTrack = genomeTrack;
            this.fileIndex = fileIndex;
            this.order = trackJson.optLong("order", 0);
        }

        boolean isMerged() {
            return type == Type.MERGED;
        }

        boolean isCombined() {
            return type == Type.COMBINED;
        }

        boolean isBlat() {
            return type == Type.BLAT;
        }

        boolean isMotif() {
            return type == Type.MOTIF;
        }

        boolean isGenome() {
            return type == Type.GENOME;
        }
    }
}
