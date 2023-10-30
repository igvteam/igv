package org.broad.igv.ucsc;

/**
 * A static json-like object, emulates javascript equivalent.   Created to ease port of session code from javascript.
 * Purpose is to hold state information,  we do not bloat this class with needless set/get methods
 */

public class TrackConfig {

    public String id;
    public String name;
    public String url;
    public String indexURL;
    public String format;
    public String displayMode;
    public String description;
    public boolean autoscale;
    public int maxHeight;
    public int height;
    public int minHeight;
    public String color;
    public String altColor;
    public int min;
    public int max;
    public boolean visible;
    public String infoURL;
    public String searchIndex;
    public String searchTrix;
    public String group;
    public int order;

    /**
     * The only required property of a track configuration is a URL (which can be an actual URL or a static file path)
     * @param url
     */
    public TrackConfig(String url) {
        this.url = url;
    }
}
