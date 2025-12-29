package org.igv.sam.mods;

import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class BaseModficationFilter {

    String modification;
    char base = 0;

    public BaseModficationFilter(String modification, char base) {
        this.modification = modification;
        this.base = base;
    }

    public BaseModficationFilter(String modification) {
        this.modification = modification;
        this.base = 0;
    }

    public boolean pass(String mod) {
        if(mod.startsWith("NONE_")) {
            char b = mod.charAt(5);
            return (this.base > 0 && this.base == b) ||
                    (modToBaseMap.containsKey(this.modification.substring(0, 1)) &&
                            modToBaseMap.get(this.modification.substring(0, 1)).charAt(0) == b);
        } else {
            return (modification == null || modification.contains(mod));
        }
    }

    public boolean pass(String mod, char b) {

        if(this.base != 0 && this.base == b) {
            return true;
        } else {
            return pass(mod);
        }
    }

    @Override
    public String toString() {
        return (modification == null ? "" : modification) + "," + (base == 0 ? "" : base);
    }

    public static BaseModficationFilter fromString(String str) {
        int idx = str.indexOf(",");
        String mod;
        char base;
        if (idx < 0) {
            // Backward compatibility
            mod = str;
            base = 0;
        } else if (idx == 0) {
            mod = null;
            base = str.charAt(1);
        } else {
            mod = str.substring(0, idx);
            base = idx == str.length() - 1 ? 0 : str.charAt(idx + 1);
        }
        return new BaseModficationFilter(mod, base);

    }

    public static Map<String, String> modToBaseMap;
    static {
        modToBaseMap = Stream.of(new String[][]{
                {"m", "C"},
                {"h", "C"},
                {"f", "C"},
                {"c", "C"},
                {"C", "C"},
                {"g", "T"},
                {"e", "T"},
                {"b", "T"},
                {"T", "T"},
                {"a", "A"},
                {"A", "A"},
                {"o", "G"},
                {"G", "G"},
                {"n", "N"},
                {"N", "N"}
        }).collect(Collectors.toMap(data -> data[0], data -> data[1]));

    }
}
