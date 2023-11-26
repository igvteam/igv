package org.broad.igv.feature.genome.load;

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
    public Boolean autoscale;
    public Integer maxHeight;
    public Integer height;
    public Integer minHeight;
    public String color;
    public String altColor;
    public Float min;
    public Float max;
    public Boolean visible;
    public String infoURL;
    public String searchIndex;
    public String searchTrix;
    public String group;
    public Integer order;
    public Integer visibilityWindow;
    public Boolean indexed;
    public Boolean hidden;

    /**
     * The only required property of a track configuration is a URL (which can be an actual URL or a static file path)
     * @param url
     */
    public TrackConfig(String url) {
        this.url = url;
    }

}
