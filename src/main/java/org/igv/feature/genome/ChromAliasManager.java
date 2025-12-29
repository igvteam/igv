package org.igv.feature.genome;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;

import java.io.IOException;
import java.util.*;

/**
 * A data/feature source helper class for managing chromosome aliasing.  Maps reference sequence names to aliases
 * used by the feature source (e.g. chr20 -> 20).
 */
public class ChromAliasManager {

    private static Logger log = LogManager.getLogger(ChromAliasManager.class);

    private Map<String, String> sequenceNames;

    private Genome genome;

    private Map<String, String> aliasCache;

    /**
     * @param sequenceNames - Sequence names defined by the data source (e.g. bam or feature file)
     * @param genome        - Reference genome object.
     */
    public ChromAliasManager(Collection<String> sequenceNames, Genome genome) {
        this.sequenceNames = new HashMap<>();
        for (String name : sequenceNames) {
            this.sequenceNames.put(name.toLowerCase(), name);
        }
        this.aliasCache = new HashMap<>();
        this.genome = genome;
    }

    /**
     * Get the alias used by the owning feature source for the given chromosome name.
     * @param chr
     * @return
     */
    public String getAliasName(String chr) {
        if (genome == null) {
            return chr;   // A no-op manager, used in testing.
        }
        try {
            if (!aliasCache.containsKey(chr)) {
                ChromAlias aliasRecord = genome.getAliasRecord(chr);
                if (aliasRecord == null) {
                    aliasCache.put(chr, null);  // No know alias, record to prevent searching again
                } else {
                    for (String alias : aliasRecord.values()) {
                        final String lowerCase = alias.toLowerCase();
                        if (sequenceNames.containsKey(lowerCase)) {
                            aliasCache.put(chr, sequenceNames.get(lowerCase));
                        }
                    }
                }
            }
            String alias = aliasCache.get(chr);
            return alias != null ? alias : chr;
        } catch (IOException e) {
            log.error("Error loading alias file ", e);
            return chr;
        }
    }

}
