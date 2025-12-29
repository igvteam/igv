package org.igv.jbrowse;

class JsonUtils {

    public static String toJson(String name, String value) {
        return "\"" + name + "\": \"" + value + "\"";
    }

    public static String toJson(String name, int value) {
        return "\"" + name + "\": " + value;
    }
}
