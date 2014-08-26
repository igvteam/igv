package org.broad.igv.ga4gh;

/**
 * Created by jrobinso on 8/25/14.
 */
public class Ga4ghProvider {

    String name;
    String baseURL;
    String authKey;

    public Ga4ghProvider(String name, String baseURL, String authKey) {
        this.name = name;
        this.baseURL = baseURL;
        this.authKey = authKey;
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
}
