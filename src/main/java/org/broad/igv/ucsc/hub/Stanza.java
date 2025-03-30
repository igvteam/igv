package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;

import java.util.*;

class Stanza {

    private static Set<String> parentOverrideProperties = new HashSet<>(Arrays.asList("visibility", "priority", "group"));
    private static Set<String> inheritableProperties = new HashSet<>(Arrays.asList("group", "priority", "color",
            "altColor", "autoscale", "viewLimits", "negativeValues", "maxWindowToQuery", "transformFun",
            "windowingFunction", "yLineMark", "yLineOnOff", "graphTypeDefault", "interactUp", "interactMultiRegion",
            "endsVisible", "maxHeightPixels", "scoreMin", "scoreFilter", "scoreFilterLimits",
            "minAliQual", "bamColorTag", "bamColorMode", "bamGrayMode", "colorByStrand", "itemRgb", "html"));
    final String type;
    final String name;
    Stanza parent;
    Map<String, String> properties;


    Stanza(String type, String name) {
        this.type = type;
        this.name = name;
        this.properties = new HashMap<>();
    }

    static String firstWord(String str) {
        return Globals.whitespacePattern.split(str)[0];
    }

    public String getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    String getOwnProperty(String key) {
        return this.properties.get(key);
    }

    boolean hasOwnProperty(String key) {
        return getProperty(key) != null;
    }

    String getProperty(String key) {

        if (this.properties.containsKey("noInherit")) {
            return this.properties.get(key);
        } else if (parentOverrideProperties.contains(key) && this.parent != null && this.parent.hasProperty(key)) {
            return this.parent.getProperty(key);
        } else if (this.properties.containsKey(key)) {
            return this.properties.get(key);
        } else if (this.parent != null && this.inheritableProperties.contains(key)) {
            return this.parent.getProperty(key);
        } else {
            return null;
        }
    }

    boolean hasProperty(String key) {
        return getProperty(key) != null;
    }

    String format() {
        String type = this.getOwnProperty("type");
        if (type != null) {
            // Trim extra bed qualifiers (e.g. bigBed + 4)
            return firstWord(type);
        }
        return null;  // unknown type
    }

    public Stanza getParent() {
        return parent;
    }

    public void setParent(Stanza parent) {
        this.parent = parent;
    }

    public Stanza getAncestor() {
        if (parent != null) {
            return parent.getAncestor();
        } else {
            return this;
        }
    }


}
