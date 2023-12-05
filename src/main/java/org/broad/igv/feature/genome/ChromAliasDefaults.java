package org.broad.igv.feature.genome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ChromAliasDefaults extends ChromAliasSource {

    public ChromAliasDefaults(String id, List<String> chromosomeNames) {

        this.init(id, chromosomeNames);
    }

    private void init(String id, List<String> chromosomeNames) {

        List<ChromAlias> aliasRecords = new ArrayList<>();
        for (String name : chromosomeNames) {

            boolean skipRest = false;
            ChromAlias record = new ChromAlias(name);
            aliasRecords.add(record);

            if (name.startsWith("gi|")) {
                // NCBI
                String alias = ChromAliasDefaults.getNCBIName(name);
                record.put("ncbi-gi-versioned", alias);

                // Also strip version number out, if present
                int dotIndex = alias.lastIndexOf('.');
                if (dotIndex > 0) {
                    alias = alias.substring(0, dotIndex);
                    record.put("ncbi-gi", alias);
                }
            } else {
                // Special cases for human and mouse
                if (id.startsWith("hg") || id.equals("1kg_ref") || id.equals("b37")) {
                    switch (name) {
                        case "chrX":
                            record.put("ncbi", "23");
                            record.put("?", "X");
                            skipRest = true;
                            break;
                        case "chrY":
                            record.put("ncbi", "24");
                            record.put("?", "Y");
                            skipRest = true;
                            break;
                    }
                } else if (id.startsWith("GRCh")) {
                    switch (name) {
                        case "23":
                            record.put("ucsc", "chrX");
                            skipRest = true;
                            break;
                        case "24":
                            record.put("ucsc", "chrY");
                            skipRest = true;
                            break;
                    }
                } else if (id.startsWith("mm") || id.startsWith("rheMac")) {
                    switch (name) {
                        case "chrX":
                            record.put("ncbi", "21");
                            skipRest = true;
                            break;
                        case "chrY":
                            record.put("ncbi", "22");
                            skipRest = true;
                            break;
                    }
                } else if (id.startsWith("GRCm")) {
                    switch (name) {
                        case "21":
                            record.put("ucsc", "chrX");
                            skipRest = true;
                            break;
                        case "22":
                            record.put("ucsc", "chrY");
                            skipRest = true;
                            break;
                    }
                }
                if (skipRest) continue;

                //
                if (name.equals("chrM")) {
                    record.put("ncbi", "MT");
                } else if (name.equals("MT")) {
                    record.put("ucsc", "chrM");
                } else if (name.toLowerCase().startsWith("chr")) {
                    record.put("ncbi", name.substring(3));
                } else if (isSmallPositiveInteger(name)) {
                    record.put("ucsc", "chr" + name);
                }
            }
        }

        for (ChromAlias rec : aliasRecords) {
            for (String a : rec.values()) {
                this.aliasCache.put(a, rec);
            }
        }
    }

    /**
     * Return the canonical chromosome name for the alias.  If none found return the alias
     *
     * @param alias
     * @returns {*}
     */
    @Override
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
    @Override
    public String getChromosomeAlias(String chr, String nameSet) {
        ChromAlias aliasRecord = this.aliasCache.get(chr);
        return aliasRecord != null && aliasRecord.containsKey(nameSet) ? aliasRecord.get(nameSet) : chr;
    }

    @Override
    public ChromAlias search(String alias) {
        return this.aliasCache.get(alias);
    }

    /**
     * Extract the user friendly name from an NCBI accession
     * example: gi|125745044|ref|NC_002229.3|  =>  NC_002229.3
     */

    public static String getNCBIName(String name) {

        String[] tokens = name.split("\\|");
        return tokens[tokens.length - 1];
    }


    public static boolean isSmallPositiveInteger(String str) {
        int length = str.length();
        if (length > 100) return false;
        for (int i = 0; i < length; i++) {
            char c = str.charAt(i);
            if (c < '0' || c > '9') {
                return false;
            }
        }
        return true;
    }
}
