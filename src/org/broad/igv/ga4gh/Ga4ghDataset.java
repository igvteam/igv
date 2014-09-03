package org.broad.igv.ga4gh;

import java.util.List;

/**
 * Created by jrobinso on 9/3/14.
 */
public class Ga4ghDataset {

    String id;
    String name;
    String genomeId;
    List<Ga4ghReadset> readsets;

    public Ga4ghDataset(String id, String name, String genomeId) {
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

    public List<Ga4ghReadset> getReadsets() {
        return readsets;
    }

    public void setReadsets(List<Ga4ghReadset> readsets) {
        this.readsets = readsets;
    }
}
