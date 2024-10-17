package org.broad.igv.feature.genome.load;

import com.google.gson.Gson;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.feature.genome.Sequence;

import java.util.LinkedHashMap;
import java.util.List;

/**
 * A static json-like object representing a genome configuration. Emulates the javascript equivalent.
 * This class was created to ease port of code from javascript.
 * <p>
 * Objects of this class are created by two paths:
 *   (1) GSON deserializtion of a genome json file, and
 *   (2) parsing track hub files.
 */

public class GenomeConfig {

    private String id;
    private String name;
    private String fastaURL;
    private String indexURL;
    private String gziIndexURL;
    private String compressedIndexURL;
    private String twoBitURL;
    private String nameSet;
    private boolean wholeGenomeView = true;
    private String defaultPos;
    private String description;
    private String blat;
    private String chromAliasBbURL;
    private String twoBitBptURL;
    private String infoURL;
    private String cytobandURL;
    private String cytobandBbURL;

    private Boolean ordered;
    private String blatDB;
    private String ucsdID;
    private String aliasURL;
    private String[] chromosomeOrder;
    private String chains;

    private String chromSizesURL;

    private List<TrackConfig> tracks;

    private List<TrackConfig> annotations;  // Backward compatibility, synonym for tracks

    // The properties below support the legacy ".genome" file, which directly loads resources from a zip archive.

    private Sequence sequence;
    private LinkedHashMap<String, List<Cytoband>> cytobands;
    private List<List<String>> chromAliases;

    public static GenomeConfig fromJson(String json) {
        return (new Gson()).fromJson(json, GenomeConfig.class);
    }

    public GenomeConfig() {
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
        return gziIndexURL;
    }

    public void setGziIndexURL(String gziIndexURL) {
        this.gziIndexURL = gziIndexURL;
    }

    public String getCompressedIndexURL() {
        return compressedIndexURL;
    }

    public void setCompressedIndexURL(String compressedIndexURL) {
        this.compressedIndexURL = compressedIndexURL;
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

    public String getUcsdID() {
        return ucsdID;
    }

    public void setUcsdID(String ucsdID) {
        this.ucsdID = ucsdID;
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

    public List<TrackConfig> getTracks() {
        return tracks;
    }

    public void setTracks(List<TrackConfig> tracks) {
        this.tracks = tracks;
    }

    public List<TrackConfig> getAnnotations() {
        return annotations;
    }

    public void setAnnotations(List<TrackConfig> annotations) {
        this.annotations = annotations;
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
}
