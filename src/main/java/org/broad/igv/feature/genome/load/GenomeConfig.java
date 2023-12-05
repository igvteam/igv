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
 * Objects of this class are created by two paths:  (1) GSON deserializtion of a genome json file, and (2)
 * parsing track hub files.
 */

public class GenomeConfig {

    public String id;
    public String name;
    public String fastaURL;
    public String indexURL;
    public String gziIndexURL;
    public String compressedIndexURL;
    public String twoBitURL;
    public String nameSet;
    public Boolean wholeGenomeView;
    public String defaultPos;
    public String description;
    public String blat;
    public String chromAliasBbURL;
    public String twoBitBptURL;
    public String infoURL;
    public String cytobandURL;
    public String cytobandBbURL;

    public Boolean ordered;
    public String blatDB;
    public String ucsdID;
    public String aliasURL;
    public String[] chromosomeOrder;
    public String chains;

    public String chromSizesURL;

    public List<TrackConfig> tracks;

    public List<TrackConfig> annotations;  // Backward compatibility, synonym for tracks

    // The properties below support the legacy ".genome" file, which directly loads resources from a zip archive.

    public Sequence sequence;
    public LinkedHashMap<String, List<Cytoband>> cytobands;
    public List<List<String>> chromAliases;

    public static GenomeConfig fromJson(String json) {
        return (new Gson()).fromJson(json, GenomeConfig.class);

    }
}
