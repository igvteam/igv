package org.broad.igv.feature.genome;

import java.io.IOException;
import java.util.*;

public class ChromAliasManager {

    Set<String> sequenceNames;
    
    Genome genome;
    
    Map<String, String> aliasCache;

    public ChromAliasManager(Collection<String> sequenceNames, Genome genome) {
        this.sequenceNames = new HashSet<>(sequenceNames);
        this.aliasCache = new HashMap<>();
        this.genome = genome;
    }

    public String getAliasName(String chr) {
        if(genome == null) {
            return chr;   // A no-op manager, used mostly in testing.
        }
        try {
            if(!aliasCache.containsKey(chr)) {
                ChromAlias aliasRecord = genome.getAliasRecord(chr);
                for(String alias : aliasRecord.values()) {
                    if(sequenceNames.contains(alias));
                    aliasCache.put(chr, alias);
                    return alias;
                }
            }
            return aliasCache.get(chr);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    
}
