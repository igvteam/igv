package org.broad.igv.feature.genome.load;

import java.util.Map;

/**
 * A static json-like object, emulates javascript equivalent.   Created to ease port of session code from javascript.
 */

public class TrackConfig implements Cloneable {

    private String id;
    private String name;
    private String longLabel;
    private String url;
    private String indexURL;
    private String trixURL;
    private String format;
    private String displayMode;
    private String description;
    private Boolean autoscale;
    private Integer maxHeight;
    private Integer height;
    private Integer minHeight;
    private String color;
    private String altColor;
    private Float min;
    private Float max;
    private Boolean visible;
    private String infoURL;
    private String searchIndex;
    private String group;
    private Integer visibilityWindow;
    private Boolean indexed;
    private Boolean hidden;
    private String html;
    private String panelName;
    private String labelField;

    private String stanzaParent;   // For supporting track hubs
    private Map<String, String> attributes;

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

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getLabelField() {
        return labelField;
    }

    public void setLabelField(String labelField) {
        this.labelField = labelField;
    }

    public String getLongLabel() {
        return longLabel;
    }

    public void setLongLabel(String longLabel) {
        this.longLabel = longLabel;
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
    }

    public String getIndexURL() {
        return indexURL;
    }

    public void setIndexURL(String indexURL) {
        this.indexURL = indexURL;
    }

    public String getFormat() {
        return format;
    }

    public void setFormat(String format) {
        this.format = format;
    }

    public String getDisplayMode() {
        return displayMode;
    }

    public void setDisplayMode(String displayMode) {
        this.displayMode = displayMode;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Boolean getAutoscale() {
        return autoscale;
    }

    public void setAutoscale(Boolean autoscale) {
        this.autoscale = autoscale;
    }

    public Integer getMaxHeight() {
        return maxHeight;
    }

    public void setMaxHeight(Integer maxHeight) {
        this.maxHeight = maxHeight;
    }

    public Integer getHeight() {
        return height;
    }

    public void setHeight(Integer height) {
        this.height = height;
    }

    public Integer getMinHeight() {
        return minHeight;
    }

    public void setMinHeight(Integer minHeight) {
        this.minHeight = minHeight;
    }

    public String getColor() {
        return color;
    }

    public void setColor(String color) {
        this.color = color;
    }

    public String getAltColor() {
        return altColor;
    }

    public void setAltColor(String altColor) {
        this.altColor = altColor;
    }

    public Float getMin() {
        return min;
    }

    public void setMin(Float min) {
        this.min = min;
    }

    public Float getMax() {
        return max;
    }

    public void setMax(Float max) {
        this.max = max;
    }

    public Boolean getVisible() {
        return visible;
    }

    public void setVisible(Boolean visible) {
        this.visible = visible;
    }

    public String getInfoURL() {
        return infoURL;
    }

    public void setInfoURL(String infoURL) {
        this.infoURL = infoURL;
    }

    public String getSearchIndex() {
        return searchIndex;
    }

    public void setSearchIndex(String searchIndex) {
        this.searchIndex = searchIndex;
    }

    public String getTrixURL() {
        return trixURL;
    }

    public void setTrixURL(String trixURL) {
        this.trixURL = trixURL;
    }

    public String getGroup() {
        return group;
    }

    public void setGroup(String group) {
        this.group = group;
    }

    public Integer getVisibilityWindow() {
        return visibilityWindow;
    }

    public void setVisibilityWindow(Integer visibilityWindow) {
        this.visibilityWindow = visibilityWindow;
    }

    public Boolean getIndexed() {
        return indexed;
    }

    public void setIndexed(Boolean indexed) {
        this.indexed = indexed;
    }

    public Boolean getHidden() {
        return hidden;
    }

    public void setHidden(Boolean hidden) {
        this.hidden = hidden;
    }

    public String getHtml() {
        return html;
    }

    public void setHtml(String html) {
        this.html = html;
    }

    public String getPanelName() {
        return panelName;
    }

    public void setPanelName(String panelName) {
        this.panelName = panelName;
    }

    public String getStanzaParent() {
        return stanzaParent;
    }

    public void setStanzaParent(String stanzaParent) {
        this.stanzaParent = stanzaParent;
    }

    @Override
    protected TrackConfig clone() throws CloneNotSupportedException {
        return (TrackConfig) super.clone();
    }

    public void setAttributes(Map<String, String> attributes) {
        this.attributes = attributes;
    }

    public Map<String, String> getAttributes() {
        return attributes;
    }

    public void setAttribute(String group, String label) {
        this.attributes.put(group, label);
    }
}
