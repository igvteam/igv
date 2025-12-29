package org.igv.feature.genome.load;

import org.igv.feature.Cytoband;
import org.igv.feature.genome.Sequence;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.json.JSONObject;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * A static json-like object representing a genome configuration. Emulates the javascript equivalent.
 * This class was created to ease port of code from javascript.
 *
 * NOTE: The property names match names in the corresponding genome json files.  They cannot be renamed.
 * <p>
 * Objects of this class are created by two paths:
 * (1) GSON deserializtion of a genome json file, and
 * (2) parsing track hub files.
 */

public class GenomeConfig implements Cloneable {

    private static Logger log = LogManager.getLogger(GenomeConfig.class);

    public String id;
    public String fastaURL;
    public String indexURL;
    public String gziIndexURL;
    public String compressedIndexURL;   // Used by reflection from GSON.  Do not remove even though IDE indicates no usages
    public String twoBitURL;
    public String twoBitBptURL;
    public String nameSet;
    public String defaultPos;
    public String description;
    public String blat;
    public String chromAliasBbURL;
    public String chromSizesURL;
    public String infoURL;
    public String cytobandURL;
    public String cytobandBbURL;
    public Boolean ordered;
    public String blatDB;
    public String ucscID;
    public String aliasURL;
    public String accession;
    public String taxId;
    public String organism;
    public String scientificName;
    public String maneBbURL;
    public String maneTrixURL;
    public String rsdbURL;
    public String[] chromosomeOrder;

    public boolean wholeGenomeView = true;
    private String name;
    private List<TrackConfig> tracks;
    private List<String> hubs;

    // The properties below support the legacy ".genome", which directly loads resources from a zip archive, as well
    // as ".gbk" files.
    private Sequence sequence;
    private LinkedHashMap<String, List<Cytoband>> cytobands;
    private List<List<String>> chromAliases;

    public static GenomeConfig fromJson(String json) {
        if (json.contains("chromosomeOrder")) {
            json = fixChromosomeOrder(json);
        }
        GenomeConfig config = new GenomeConfig();
        JSONObject jsonObj = new JSONObject(json);
        config.id = jsonObj.optString("id");
        config.name = jsonObj.optString("name", null);
        config.fastaURL = jsonObj.optString("fastaURL", null);
        config.indexURL = jsonObj.optString("indexURL", null);
        config.gziIndexURL = jsonObj.optString("gziIndexURL", null);
        config.compressedIndexURL = jsonObj.optString("compressedIndexURL", null);
        config.twoBitURL = jsonObj.optString("twoBitURL", null);
        config.twoBitBptURL = jsonObj.optString("twoBitBptURL", null);
        config.nameSet = jsonObj.optString("nameSet", null);
        config.defaultPos = jsonObj.optString("defaultPos", null);
        config.description = jsonObj.optString("description", null);
        config.blat = jsonObj.optString("blat", null);
        config.chromAliasBbURL = jsonObj.optString("chromAliasBbURL", null);
        config.chromSizesURL = jsonObj.optString("chromSizesURL", null);
        config.infoURL = jsonObj.optString("infoURL", null);
        config.cytobandURL = jsonObj.optString("cytobandURL", null);
        config.cytobandBbURL = jsonObj.optString("cytobandBbURL", null);
        config.ordered = jsonObj.has("ordered") ? jsonObj.optBoolean("ordered") : null;
        config.blatDB = jsonObj.optString("blatDB", null);
        config.ucscID = jsonObj.optString("ucscID", null);
        config.aliasURL = jsonObj.optString("aliasURL", null);
        config.accession = jsonObj.optString("accession", null);
        config.taxId = jsonObj.optString("taxId", null);
        config.organism = jsonObj.optString("organism", null);
        config.scientificName = jsonObj.optString("scientificName", null);
        config.maneBbURL = jsonObj.optString("maneBbURL", null);
        config.maneTrixURL = jsonObj.optString("maneTrixURL", null);
        config.rsdbURL = jsonObj.optString("rsdbURL", null);
        config.chromosomeOrder = jsonObj.has("chromosomeOrder") ?
                jsonObj.getJSONArray("chromosomeOrder").toList().toArray(new String[0]) : null;

        config.wholeGenomeView = jsonObj.optBoolean("wholeGenomeView", true);

        config.hubs = new ArrayList<>();
        if (jsonObj.has("hubs")) {
            for (int i = 0; i < jsonObj.getJSONArray("hubs").length(); i++) {
                config.hubs.add(jsonObj.getJSONArray("hubs").getString(i));
            }
        }
        if (jsonObj.has("tracks")) {
            config.tracks = new ArrayList<>();
            for (int i = 0; i < jsonObj.getJSONArray("tracks").length(); i++) {
                JSONObject track = jsonObj.getJSONArray("tracks").getJSONObject(i);
                config.tracks.add(TrackConfig.fromJSON(track.toString()));
            }
        }

        return config;
    }

    public JSONObject toJSON() {
        JSONObject jsonObj = new JSONObject();

        if (id != null) jsonObj.put("id", id);
        if (name != null) jsonObj.put("name", name);
        if (fastaURL != null) jsonObj.put("fastaURL", fastaURL);
        if (indexURL != null) jsonObj.put("indexURL", indexURL);
        if (gziIndexURL != null) jsonObj.put("gziIndexURL", gziIndexURL);
        if (compressedIndexURL != null) jsonObj.put("compressedIndexURL", compressedIndexURL);
        if (twoBitURL != null) jsonObj.put("twoBitURL", twoBitURL);
        if (twoBitBptURL != null) jsonObj.put("twoBitBptURL", twoBitBptURL);
        if (nameSet != null) jsonObj.put("nameSet", nameSet);
        if (defaultPos != null) jsonObj.put("defaultPos", defaultPos);
        if (description != null) jsonObj.put("description", description);
        if (blat != null) jsonObj.put("blat", blat);
        if (chromAliasBbURL != null) jsonObj.put("chromAliasBbURL", chromAliasBbURL);
        if (chromSizesURL != null) jsonObj.put("chromSizesURL", chromSizesURL);
        if (infoURL != null) jsonObj.put("infoURL", infoURL);
        if (cytobandURL != null) jsonObj.put("cytobandURL", cytobandURL);
        if (cytobandBbURL != null) jsonObj.put("cytobandBbURL", cytobandBbURL);
        if (ordered != null) jsonObj.put("ordered", ordered);
        if (blatDB != null) jsonObj.put("blatDB", blatDB);
        if (ucscID != null) jsonObj.put("ucscID", ucscID);
        if (aliasURL != null) jsonObj.put("aliasURL", aliasURL);
        if (accession != null) jsonObj.put("accession", accession);
        if (taxId != null) jsonObj.put("taxId", taxId);
        if (organism != null) jsonObj.put("organism", organism);
        if (scientificName != null) jsonObj.put("scientificName", scientificName);
        if (chromosomeOrder != null) jsonObj.put("chromosomeOrder", chromosomeOrder);
        if (maneBbURL != null) jsonObj.put("maneBbURL", maneBbURL);
        if (maneTrixURL != null) jsonObj.put("maneTrixURL", maneTrixURL);
        if (rsdbURL != null) jsonObj.put("rsdbURL", rsdbURL);
        jsonObj.put("wholeGenomeView", wholeGenomeView);

        if (hubs != null && !hubs.isEmpty()) jsonObj.put("hubs", hubs);

        if (tracks != null && !tracks.isEmpty()) {
            List<JSONObject> trackJsons = new ArrayList<>();
            for (TrackConfig track : tracks) {
                trackJsons.add(track.toJSON());
            }
            jsonObj.put("tracks", trackJsons);
        }

        return jsonObj;
    }

    public GenomeConfig() {
    }

    public List<String> getHubs() {
        return hubs;
    }

    public void setHubs(List<String> hubs) {
        this.hubs = hubs;
    }

    public void setName(String name) {
        this.name = name;
    }

    /**
     * If 'name' is not explicitly set derive it from known properties.
     *
     * @return
     */
    public String getName() {
        if (name == null) {
            if (this.scientificName != null) {
                name = this.scientificName;
            } else if (this.organism != null) {
                name = this.organism;
            } else if (this.description != null) {
                name = this.description;
            }
            if (name == null) {
                name = id;
            } else {
                name = name + " (" + id + ")";
            }
        }
        return name;
    }

    public String getUcscID() {
        return ucscID == null ? id : ucscID;
    }

    public String[] getChromosomeOrder() {
        return chromosomeOrder;
    }

    public List<TrackConfig> getTrackConfigs() {
        return tracks;
    }

    public void setTracks(List<TrackConfig> tracks) {
        this.tracks = tracks;
    }

    public Sequence getSequence() {
        return sequence;
    }

    public void setSequence(Sequence sequence) {
        this.sequence = sequence;
    }

    public LinkedHashMap<String, List<Cytoband>> getCytobands() {
        return cytobands;
    }

    public void setCytobands(LinkedHashMap<String, List<Cytoband>> cytobands) {
        this.cytobands = cytobands;
    }

    public List<List<String>> getChromAliases() {
        return chromAliases;
    }

    public void setChromAliases(List<List<String>> chromAliases) {
        this.chromAliases = chromAliases;
    }

    public GenomeConfig copy() {
        return this.clone();
    }

    protected GenomeConfig clone() {
        try {
            GenomeConfig clone = (GenomeConfig) super.clone();

            if (this.tracks != null) {
                clone.tracks = new ArrayList<>();
                for (TrackConfig trackConfig : this.tracks) {
                    clone.tracks.add(trackConfig.clone());
                }
            }

            return clone;
        } catch (CloneNotSupportedException e) {
            log.error("Cloning not supported", e);
            return this;
        }
    }

    /**
     * Fix deprecated form of chromosome order (comma delimited list of strings)
     *
     * @param jsonString
     * @return
     */
    private static String fixChromosomeOrder(String jsonString) {
        JSONObject obj = new JSONObject(jsonString);
        Object chromosomeOrder = obj.opt("chromosomeOrder");
        if (chromosomeOrder != null) {
            if (chromosomeOrder instanceof String) {
                obj.put("chromosomeOrder", Arrays.stream(((String) chromosomeOrder).split(",")).map(String::trim).toArray(String[]::new));
            }
        }
        return obj.toString();
    }

    public void removeHub(String url) {
        this.hubs.remove(url);
    }
}
