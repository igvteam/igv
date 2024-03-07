package org.broad.igv.feature.genome;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.io.IOException;
import java.util.*;

/**
 * A data/feature source helper class for managing chromosome aliasing.  Maps reference sequence names to aliases
 * used by the feature source (e.g. chr20 -> 20).
 */
public class ChromAliasManager {

    private static Logger log = LogManager.getLogger(ChromAliasManager.class)

    private Set<String> sequenceNames;

    private Genome genome;

    private Map<String, String> aliasCache;

    /**
     * @param sequenceNames - Sequence names defined by the data source (e.g. bam or feature file)
     * @param genome        - Reference genome object.
     */
    public ChromAliasManager(Collection<String> sequenceNames, Genome genome) {
        this.sequenceNames = new HashSet<>(sequenceNames);
        this.aliasCache = new HashMap<>();
        this.genome = genome;
    }

    public String getAliasName(String seqName) {
        if (genome == null) {
            return seqName;   // A no-op manager, used in testing.
        }
        try {
            if (!aliasCache.containsKey(seqName)) {
                ChromAlias aliasRecord = genome.getAliasRecord(seqName);
                if (aliasRecord == null) {
                    aliasCache.put(seqName, null);  // No know alias, record to prevent searching again
                } else {
                    for (String alias : aliasRecord.values()) {
                        if (sequenceNames.contains(alias)) {
                            aliasCache.put(seqName, alias);
                        }
                    }
                }
            }
            String alias = aliasCache.get(seqName);
            return alias != null ? alias : seqName;
        } catch (IOException e) {
            log.error("Error loading alias file ", e);
            return seqName;
        }
    }

}
