package org.broad.igv.feature.genome;

import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Represents a flat text chromosome alias file, either legacy IGV or UCSC
 */
public class ChromAliasFile extends ChromAliasSource {

    private String[] nameSets = null;

    public ChromAliasFile(String path, List<String> chromosomeNames) throws IOException {
        try (BufferedReader br = ParsingUtils.openBufferedReader(path)) {

            this.aliasCache = new HashMap<>();

            Set<String> chromosomeNameSet = chromosomeNames != null ? new HashSet<>(chromosomeNames) : Collections.EMPTY_SET;

            // First line
            // NOTE -- different conventions have been used for the first line.  Older chr alias files do not specify
            // name sets.  So we apply some heuristics.
            //   (1) no column header contains internal spaces (multiple words)
            //   (2) # of header columns == # of columns in data rows
            //
            // Examples of non-conforming and conforming alias file headers:
            //   # sequenceName	alias names	UCSC database: mm10
            //   # refseq	assembly	genbank	ncbi	ucsc

            boolean firstLine = true;
            String [] headers = null;
            String line;
            while ((line = br.readLine()) != null) {
                if (firstLine && line.startsWith("#")) {
                    String[] tokens = line.substring(1).trim().split("\\s*\t\\s*");
                    if(Arrays.stream(tokens).noneMatch(token -> token.contains(" "))) {
                        headers = tokens;
                    }
                } else {
                    String[] tokens = line.split("\t");
                    if (tokens.length > 1) {

                        if(nameSets == null && headers != null) {
                            if(headers.length == tokens.length) {
                                nameSets = headers;
                            }
                        }

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
    }

    /**
     * Constructor for legacy  .genome files
     *
     * @throws IOException
     */
    public ChromAliasFile(List<List<String>> chromAliases, List<String> chromosomeNames) throws IOException {

        Set<String> chromosomeNameSet = chromosomeNames != null ? new HashSet<>(chromosomeNames) : Collections.EMPTY_SET;

        for (List<String> aliases : chromAliases) {
            if (aliases.size() > 1) {

                // Find the canonical chromosome
                String chr = null;
                for (String c : aliases) {
                    if (chromosomeNameSet.contains(c)) {
                        chr = c;
                        break;
                    }
                }
                if (chr == null) {
                    chr = aliases.get(0);
                }

                ChromAlias aliasRecord = new ChromAlias(chr);
                this.aliasCache.put(chr, aliasRecord);
                int idx = 0;
                for (String a : aliases) {
                    String key = String.valueOf(idx++);
                    aliasRecord.put(key, a);
                    this.aliasCache.put(a, aliasRecord);
                }
            }
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


    public ChromAlias search(String alias) {
        return this.aliasCache.get(alias);
    }

    public boolean hasNameSets() {
        return nameSets != null;
    }
}