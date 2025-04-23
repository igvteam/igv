package org.broad.igv.ucsc.hub;


public class HubDescriptor {

    private String url;
    private String shortLabel;
    private String longLabel;

    public HubDescriptor(String shortLabel, String longLabel, String url) {
        this.longLabel = longLabel;
        this.shortLabel = shortLabel;
        this.url = url;
    }

    public String getUrl() {
        return url;
    }

    public String getShortLabel() {
        return shortLabel;
    }

    public String getLongLabel() {
        return longLabel;
    }

}
