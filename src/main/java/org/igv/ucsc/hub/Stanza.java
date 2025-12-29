package org.igv.ucsc.hub;

import org.igv.Globals;

import java.util.*;

class Stanza {

    private static Set<String> parentOverrideProperties = new HashSet<>(Arrays.asList("visibility", "priority", "group"));
    private static final Set<String> nonInheritableProperties = new HashSet<>(Arrays.asList(
            "track", "type", "shortLabel", "longLabel", "bigDataUrl",
            "parent", "superTrack", "priority", "compositeTrack", "view", "compositeContainer"
    ));

    final String type;
    final String name;
    private Stanza parent;
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
            return this.getOwnProperty(key);
        } else if (parentOverrideProperties.contains(key) && this.parent != null && this.parent.hasProperty(key)) {
            return this.parent.getProperty(key);
        } else if (this.properties.containsKey(key)) {
            return this.properties.get(key);
        } else if (this.parent != null && !this.nonInheritableProperties.contains(key)) {
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
        if (parent == this) {
            throw new IllegalArgumentException("Stanza cannot be its own parent");
        }
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
