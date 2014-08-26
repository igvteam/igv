package org.broad.igv.ga4gh;

/**
 * Minimal representation of a readset for prototype purposes.
 *
 */
public class Ga4ghReadset {

    String id;
    String name;
    String genomeId;

    public Ga4ghReadset(String id, String name, String genomeId) {
        this.id = id;
        this.name = name;
        this.genomeId = genomeId;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public String getGenomeId() {
        return genomeId;
    }
}
