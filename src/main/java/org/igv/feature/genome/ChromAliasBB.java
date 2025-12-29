package org.igv.feature.genome;

import org.igv.feature.BasicFeature;
import org.igv.feature.IGVFeature;
import org.igv.ucsc.bb.BBFile;

import java.io.IOException;
import java.util.List;

public class ChromAliasBB extends ChromAliasSource {


    BBFile reader;

    public ChromAliasBB(String path, Genome genome) throws IOException {
        this.reader = new BBFile(path, genome);
    }

    void preload(List<String> chromosomeNames) {
        try {
            this.reader.preload();
            for (String chr : chromosomeNames) {
                search(chr);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Return an alternate chromosome name (alias).
     *
     * @param chr
     * @param nameSet -- The name set, e.g. "ucsc"
     * @returns {*|undefined}
     */
    public String getChromosomeAlias(String chr, String nameSet) {
        ChromAlias aliasRecord = this.aliasCache.get(chr);
        return aliasRecord != null && aliasRecord.containsKey(nameSet) ? aliasRecord.get(nameSet) : chr;
    }

    /**
     * Search for chromosome alias bed record.  If found, cache results in the alias -> chr map
     *
     * @param alias
     * @returns {Promise<any>}
     */
    public ChromAlias search(String alias) throws IOException {
        if (!this.aliasCache.containsKey(alias)) {
            List<IGVFeature> results = this.reader.search(alias);
            if (results != null) {
                for (IGVFeature f : results) {
                    String chr = f.getChr();
                    ChromAlias aliasRecord = new ChromAlias(chr);
                    this.aliasCache.put(chr, aliasRecord);
                    for (String key : f.getAttributeKeys()) {
                        final String a = f.getAttribute(key);
                        aliasRecord.put(key, a);
                        this.aliasCache.put(a, aliasRecord);      // One entry for each alias
                    }
                }
            }
        }
        return this.aliasCache.get(alias);
    }
}