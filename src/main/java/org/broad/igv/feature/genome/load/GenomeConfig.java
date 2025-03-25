package org.broad.igv.feature.genome.load;

import com.google.gson.Gson;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.util.*;

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

    private String id;
    private String name;
    private String fastaURL;
    private String indexURL;
    private String gziIndexURL;
    private String compressedIndexURL;   // Used by reflection from GSON.  Do not remove
    private String twoBitURL;
    private String twoBitBptURL;
    private String nameSet;
    private boolean wholeGenomeView = true;
    private String defaultPos;
    private String description;
    private String blat;
    private String chromAliasBbURL;

    private String infoURL;
    private String cytobandURL;
    private String cytobandBbURL;

    private Boolean ordered;
    private String blatDB;
    private String ucscID;
    private String aliasURL;
    private String[] chromosomeOrder;
    private String chains;

    private String chromSizesURL;

    private List<TrackConfig> tracks;

    private List<String> hubs;

    // The properties below support the legacy ".genome" file, which directly loads resources from a zip archive.

    private Sequence sequence;
    private LinkedHashMap<String, List<Cytoband>> cytobands;
    private List<List<String>> chromAliases;
    private boolean fromJson = false;  // Until proven otherwise

    public static GenomeConfig fromJson(String json) {
        if (json.contains("chromosomeOrder")) {
            json = fixChromosomeOrder(json);
        }
        GenomeConfig config = (new Gson()).fromJson(json, GenomeConfig.class);
        config.fromJson = true;
        return config;
    }

    public GenomeConfig() {
    }

    public String toJson() {
        return new Gson().toJson(this);
    }

    public List<String> getHubs() {
        return hubs;
    }

    public void setHubs(List<String> hubs) {
        this.hubs = hubs;
    }

    public void addHub(String hub) {
        if (hubs == null) {
            hubs = new ArrayList();
        }
        hubs.add(hub);
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getFastaURL() {
        return fastaURL;
    }

    public void setFastaURL(String fastaURL) {
        this.fastaURL = fastaURL;
    }

    public String getIndexURL() {
        return indexURL;
    }

    public void setIndexURL(String indexURL) {
        this.indexURL = indexURL;
    }

    public String getGziIndexURL() {
        return this.gziIndexURL != null ? gziIndexURL : compressedIndexURL;
    }

    public void setGziIndexURL(String gziIndexURL) {
        this.gziIndexURL = gziIndexURL;
    }

    public String getTwoBitURL() {
        return twoBitURL;
    }

    public void setTwoBitURL(String twoBitURL) {
        this.twoBitURL = twoBitURL;
    }

    public String getNameSet() {
        return nameSet;
    }

    public void setNameSet(String nameSet) {
        this.nameSet = nameSet;
    }

    public boolean isWholeGenomeView() {
        return wholeGenomeView;
    }

    public void setWholeGenomeView(boolean wholeGenomeView) {
        this.wholeGenomeView = wholeGenomeView;
    }

    public String getDefaultPos() {
        return defaultPos;
    }

    public void setDefaultPos(String defaultPos) {
        this.defaultPos = defaultPos;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getBlat() {
        return blat;
    }

    public void setBlat(String blat) {
        this.blat = blat;
    }

    public String getChromAliasBbURL() {
        return chromAliasBbURL;
    }

    public void setChromAliasBbURL(String chromAliasBbURL) {
        this.chromAliasBbURL = chromAliasBbURL;
    }

    public String getTwoBitBptURL() {
        return twoBitBptURL;
    }

    public void setTwoBitBptURL(String twoBitBptURL) {
        this.twoBitBptURL = twoBitBptURL;
    }

    public String getInfoURL() {
        return infoURL;
    }

    public void setInfoURL(String infoURL) {
        this.infoURL = infoURL;
    }

    public String getCytobandURL() {
        return cytobandURL;
    }

    public void setCytobandURL(String cytobandURL) {
        this.cytobandURL = cytobandURL;
    }

    public String getCytobandBbURL() {
        return cytobandBbURL;
    }

    public void setCytobandBbURL(String cytobandBbURL) {
        this.cytobandBbURL = cytobandBbURL;
    }

    public Boolean getOrdered() {
        return ordered;
    }

    public void setOrdered(Boolean ordered) {
        this.ordered = ordered;
    }

    public String getBlatDB() {
        return blatDB;
    }

    public void setBlatDB(String blatDB) {
        this.blatDB = blatDB;
    }

    public String getUcscID() {
        return ucscID == null ? id : ucscID;
    }
    public void setUcscID(String ucscID) {
        this.ucscID = ucscID;
    }

    public String getAliasURL() {
        return aliasURL;
    }

    public void setAliasURL(String aliasURL) {
        this.aliasURL = aliasURL;
    }

    public String[] getChromosomeOrder() {
        return chromosomeOrder;
    }

    public void setChromosomeOrder(String[] chromosomeOrder) {
        this.chromosomeOrder = chromosomeOrder;
    }

    public String getChains() {
        return chains;
    }

    public void setChains(String chains) {
        this.chains = chains;
    }

    public String getChromSizesURL() {
        return chromSizesURL;
    }

    public void setChromSizesURL(String chromSizesURL) {
        this.chromSizesURL = chromSizesURL;
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
        Map obj = (new Gson()).fromJson(jsonString, Map.class);
        Object chromosomeOrder = obj.get("chromosomeOrder");
        if (chromosomeOrder != null) {
            if (chromosomeOrder instanceof String) {
                obj.put("chromosomeOrder", Arrays.stream(((String) chromosomeOrder).split(",")).map(c -> c.trim()).toArray());
            }
        }
        return (new Gson()).toJson(obj);
    }

    public boolean isFromJson() {
        return fromJson;
    }
}
