package org.broad.igv.feature.genome;

import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Represents a flat text chromosome alias file, either legacy IGV or UCSC
 */
public class ChromAliasFile extends ChromAliasSource {

    private String[] nameSets;


    public ChromAliasFile(String path, Genome genome) throws IOException {
        try (BufferedReader br = ParsingUtils.openBufferedReader(path)) {
            init(br, genome.getAllChromosomeNames());
        }
    }

    /**
     * Constructor for legacy  .genome files
     * @param br
     * @param genome
     * @throws IOException
     */
    public ChromAliasFile(BufferedReader br, Genome genome) throws IOException {
            init(br, genome.getAllChromosomeNames());

    }

    private void init(BufferedReader br, List<String> chromosomeNames) throws IOException {

        this.aliasCache = new HashMap<>();

        Set<String> chromosomeNameSet = chromosomeNames != null ? new HashSet<>(chromosomeNames) : Collections.EMPTY_SET;

        // First line
        boolean firstLine = true;
        String line;
        while ((line = br.readLine()) != null) {
            if (firstLine && line.startsWith("#")) {
                String[] tokens = line.substring(1).split("\t");
                this.nameSets = new String[tokens.length];
                for (int i = 0; i < tokens.length; i++) {
                    this.nameSets[i] = tokens[i].trim();
                }
            } else {
                String[] tokens = line.split("\t");
                if (tokens.length > 1) {

                    // Find the canonical chromosome
                    String chr = null;
                    for (String c : tokens) {
                        if (chromosomeNameSet.contains(c)) {
                            chr = c;
                            break;
                        }
                    }
                    if (chr == null) {
                        chr = tokens[0];
                    }

                    ChromAlias aliasRecord = new ChromAlias(chr);
                    this.aliasCache.put(chr, aliasRecord);
                    for (int i = 0; i < tokens.length; i++) {
                        String key = this.nameSets != null ? this.nameSets[i] : String.valueOf(i);
                        aliasRecord.put(key, tokens[i]);
                        this.aliasCache.put(tokens[i], aliasRecord);
                    }
                }
            }
            firstLine = false;
        }

    }

    /**
     * Return the canonical chromosome name for the alias.  If none found return the alias
     *
     * @param alias
     * @returns {*}
     */
    public String getChromosomeName(String alias) {
        return this.aliasCache.containsKey(alias) ? this.aliasCache.get(alias).get("chr") : alias;
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


    public ChromAlias search(String alias) {
        return this.aliasCache.get(alias);
    }

}