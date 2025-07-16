package org.broad.igv.feature.genome.load;

import java.util.Map;

/**
 * A static json-like object, emulates javascript equivalent.   Created to ease port of session code from javascript.
 */

public class TrackConfig implements Cloneable {

    public String id;
    public String name;
    public String longLabel;
    public String url;
    public String indexURL;
    public String trixURL;
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
    public String group;
    public Integer visibilityWindow;
    public Boolean indexed;
    public Boolean hidden;
    public String html;
    public String panelName;
    public String labelField;
    public String autoscaleGroup;
    public String [] filterTypes;
    public String stanzaParent; // For supporting track hubs
    public Map<String, String> attributes;

    public TrackConfig() {
    }

    /**
     * The only required property of a track configuration is a URL (which can be an actual URL or a static file path)
     *
     * @param url
     */
    public TrackConfig(String url) {
        this.url = url;
    }

    @Override
    protected TrackConfig clone() throws CloneNotSupportedException {
        return (TrackConfig) super.clone();
    }


    public void setAttribute(String group, String label) {
        this.attributes.put(group, label);
    }
}
