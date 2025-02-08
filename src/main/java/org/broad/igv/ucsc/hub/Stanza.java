package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;

import java.util.*;

class Stanza {

    private static Set<String> parentOverrideProperties = new HashSet<>(Arrays.asList("visibility", "priority", "group"));
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

    void setProperty(String key, String value) {
        this.properties.put(key, value);
    }

    String getOwnProperty(String key) {
        return this.properties.get(key);
    }

    String getProperty(String key) {

        if (parentOverrideProperties.contains(key) && this.parent != null && this.parent.hasProperty(key)) {
            return this.parent.getProperty(key);
        } else if (this.properties.containsKey(key)) {
            return this.properties.get(key);
        } else if (this.parent != null) {
            return this.parent.getProperty(key);
        } else {
            return null;
        }
    }

    boolean hasProperty(String key) {
        if (this.properties.containsKey(key)) {
            return true;
        } else if (this.parent != null) {
            return this.parent.hasProperty(key);
        } else {
            return false;
        }
    }

    String format() {
        String type = this.getProperty("type");
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
}
