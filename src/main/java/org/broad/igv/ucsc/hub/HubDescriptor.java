package org.broad.igv.ucsc.hub;


public class HubDescriptor {

    private String url;
    private String shortLabel;
    private String longLabel;
    private String descriptionUrl;
    private String dbList;

    public HubDescriptor(String shortLabel, String longLabel, String url, String descriptionUrl, String dbList) {
        this.longLabel = longLabel;
        this.shortLabel = shortLabel;
        this.url = url;
        this.descriptionUrl = descriptionUrl;
        this.dbList = dbList;
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

    public String getDescriptionUrl() {
        return descriptionUrl;
    }

    public String getDbList() {
        return dbList;
    }

    public String toString() {
        return "HubDescriptor{" +
                "url='" + url + '\'' +
                ", shortLabel='" + shortLabel + '\'' +
                ", longLabel='" + longLabel + '\'' +
                ", descriptionUrl='" + descriptionUrl + '\'' +
                ", dbList='" + dbList + '\'' +
                '}';
    }

}
