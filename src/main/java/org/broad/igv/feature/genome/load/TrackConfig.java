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
    public String[] filterTypes;
    public String stanzaParent; // For supporting track hubs
    public Map<String, String> attributes;

    public TrackConfig() {
    }

    /**
     * Create a TrackConfig from a JSON string.
     *
     * @param json JSON string representing a track configuration.
     * @return A TrackConfig object populated with the data from the JSON string.
     */
    public static TrackConfig fromJSON(String json) {

        org.json.JSONObject jsonObject = new org.json.JSONObject(json);
        TrackConfig config = new TrackConfig();
        config.id = jsonObject.optString("id", null);
        config.name = jsonObject.optString("name", null);
        config.longLabel = jsonObject.optString("longLabel", null);
        config.url = jsonObject.optString("url", null);
        config.indexURL = jsonObject.optString("indexURL", null);
        config.trixURL = jsonObject.optString("trixURL", null);
        config.format = jsonObject.optString("format", null);
        config.displayMode = jsonObject.optString("displayMode", null);
        config.description = jsonObject.optString("description", null);
        config.autoscale = jsonObject.has("autoscale") ? jsonObject.optBoolean("autoscale") : null;
        config.maxHeight = jsonObject.has("maxHeight") ? jsonObject.optInt("maxHeight") : null;
        config.height = jsonObject.has("height") ? jsonObject.optInt("height") : null;
        config.minHeight = jsonObject.has("minHeight") ? jsonObject.optInt("minHeight") : null;
        config.color = jsonObject.optString("color", null);
        config.altColor = jsonObject.optString("altColor", null);
        config.min = jsonObject.has("min") ? (float) jsonObject.optDouble("min") : null;
        config.max = jsonObject.has("max") ? (float) jsonObject.optDouble("max") : null;
        config.visible = jsonObject.has("visible") ? jsonObject.optBoolean("visible") : null;
        config.infoURL = jsonObject.optString("infoURL", null);
        config.searchIndex = jsonObject.optString("searchIndex", null);
        config.group = jsonObject.optString("group", null);
        config.visibilityWindow = jsonObject.has("visibilityWindow") ? jsonObject.optInt("visibilityWindow") : null;
        config.indexed = jsonObject.has("indexed") ? jsonObject.optBoolean("indexed") : null;
        config.hidden = jsonObject.has("hidden") ? jsonObject.optBoolean("hidden") : null;
        config.html = jsonObject.optString("html", null);
        config.panelName = jsonObject.optString("panelName", null);
        config.labelField = jsonObject.optString("labelField", null);
        config.autoscaleGroup = jsonObject.optString("autoscaleGroup", null);
        if (jsonObject.has("filterTypes")) {
            org.json.JSONArray arr = jsonObject.optJSONArray("filterTypes");
            if (arr != null) {
                config.filterTypes = new String[arr.length()];
                for (int i = 0; i < arr.length(); i++) {
                    config.filterTypes[i] = arr.optString(i, null);
                }
            }
        }
        config.stanzaParent = jsonObject.optString("stanzaParent", null);
        if (jsonObject.has("attributes")) {
            org.json.JSONObject attrs = jsonObject.optJSONObject("attributes");
            if (attrs != null) {
                config.attributes = new java.util.HashMap<>();
                for (String key : attrs.keySet()) {
                    config.attributes.put(key, attrs.optString(key, null));
                }
            }
        }
        return config;
    }


    /**
     * Convert the TrackConfig object to a JSON string.
     * <p>
     * @return A JSON string representing the TrackConfig object.
     */
    public org.json.JSONObject toJSON() {
        org.json.JSONObject jsonObject = new org.json.JSONObject();
        if (id != null) jsonObject.put("id", id);
        if (name != null) jsonObject.put("name", name);
        if (longLabel != null) jsonObject.put("longLabel", longLabel);
        if (url != null) jsonObject.put("url", url);
        if (indexURL != null) jsonObject.put("indexURL", indexURL);
        if (trixURL != null) jsonObject.put("trixURL", trixURL);
        if (format != null) jsonObject.put("format", format);
        if (displayMode != null) jsonObject.put("displayMode", displayMode);
        if (description != null) jsonObject.put("description", description);
        if (autoscale != null) jsonObject.put("autoscale", autoscale);
        if (maxHeight != null) jsonObject.put("maxHeight", maxHeight);
        if (height != null) jsonObject.put("height", height);
        if (minHeight != null) jsonObject.put("minHeight", minHeight);
        if (color != null) jsonObject.put("color", color);
        if (altColor != null) jsonObject.put("altColor", altColor);
        if (min != null) jsonObject.put("min", min);
        if (max != null) jsonObject.put("max", max);
        if (visible != null) jsonObject.put("visible", visible);
        if (infoURL != null) jsonObject.put("infoURL", infoURL);
        if (searchIndex != null) jsonObject.put("searchIndex", searchIndex);
        if (group != null) jsonObject.put("group", group);
        if (visibilityWindow != null) jsonObject.put("visibilityWindow", visibilityWindow);
        if (indexed != null) jsonObject.put("indexed", indexed);
        if (hidden != null) jsonObject.put("hidden", hidden);
        if (html != null) jsonObject.put("html", html);
        if (panelName != null) jsonObject.put("panelName", panelName);
        if (labelField != null) jsonObject.put("labelField", labelField);
        if (autoscaleGroup != null) jsonObject.put("autoscaleGroup", autoscaleGroup);
        if (filterTypes != null) {
            org.json.JSONArray arr = new org.json.JSONArray();
            for (String type : filterTypes) {
                if (type != null) arr.put(type);
            }
            jsonObject.put("filterTypes", arr);
        }
        if (stanzaParent != null) jsonObject.put("stanzaParent", stanzaParent);
        if (attributes != null) {
            org.json.JSONObject attrs = new org.json.JSONObject();
            for (Map.Entry<String, String> entry : attributes.entrySet()) {
                if (entry.getKey() != null && entry.getValue() != null) {
                    attrs.put(entry.getKey(), entry.getValue());
                }
            }
            jsonObject.put("attributes", attrs);
        }
        return jsonObject;
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
