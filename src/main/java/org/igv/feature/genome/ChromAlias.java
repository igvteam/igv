package org.igv.feature.genome;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class ChromAlias {

    private String chr;

    /**
     * Map of name set -> alias
     */
    private Map<String, String> aliases;

    public ChromAlias(String chr) {
        this.chr = chr;
        this.aliases = new HashMap<>();
        String defaultAlias = chr.startsWith("chr") ? chr.substring(3) : "chr" + chr;
        aliases.put("**DEFAULT**", defaultAlias);

    }

    public String getChr() {
        return chr;
    }

    public void put(String nameSet, String alias) {
        aliases.put(nameSet, alias);
    }

    public String get(String nameSet) {
        return aliases.get(nameSet);
    }

    public boolean containsKey(String nameSet) {
        return aliases.containsKey(nameSet);
    }

    public Collection<String> values() {
        return aliases.values();
    }


}
