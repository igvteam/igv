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
        GenomeConfig config = (new Gson()).fromJson(json, GenomeConfig.class);
        return config;
    }

    public GenomeConfig() {
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
        Map obj = (new Gson()).fromJson(jsonString, Map.class);
        Object chromosomeOrder = obj.get("chromosomeOrder");
        if (chromosomeOrder != null) {
            if (chromosomeOrder instanceof String) {
                obj.put("chromosomeOrder", Arrays.stream(((String) chromosomeOrder).split(",")).map(c -> c.trim()).toArray());
            }
        }
        return (new Gson()).toJson(obj);
    }

    public void removeHub(String url) {
        this.hubs.remove(url);
    }
}
