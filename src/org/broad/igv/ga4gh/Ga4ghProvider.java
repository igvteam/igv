package org.broad.igv.ga4gh;

import java.util.List;

/**
 * Created by jrobinso on 8/25/14.
 */
public class Ga4ghProvider {

    String name;
    String baseURL;
    String authKey;
    List<Ga4ghDataset> datasets;
    private Object id;

    public Ga4ghProvider(String name, String baseURL, String authKey, List<Ga4ghDataset> datasets) {
        this.name = name;
        this.baseURL = baseURL;
        this.authKey = authKey;
        this.datasets = datasets;
    }

    public String getName() {
        return name;
    }

    public String getBaseURL() {
        return baseURL;
    }

    public String getAuthKey() {
        return authKey;
    }

    public List<Ga4ghDataset> getDatasets() {
        return datasets;
    }

    public String toString() {
        return name;
    }

    public Object getId() {
        return id;
    }
}
