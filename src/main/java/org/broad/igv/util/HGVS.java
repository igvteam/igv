// Java
package org.broad.igv.util;

import org.broad.igv.feature.genome.ChromAlias;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.action.SearchCommand;

import java.io.IOException;

public class HGVS {

    /**
     * Returns genomic HGVS notation: <RefSeqAccession>:g.<position>
     * Example: NC_000001.11:g.1234567
     */
    public static String getHGVSNotation(Genome genome, String chr, int position) {
        try {
            ChromAlias aliasRecord = genome.getAliasRecord(chr);
            String accession = null;

            if (aliasRecord != null) {
                for (String alias : aliasRecord.values()) {
                    if (alias.startsWith("NC_")) {
                        accession = alias;
                        break;
                    }
                }
            }

            // Fallback to provided chromosome if no RefSeq accession is found
            if (accession == null || accession.isEmpty()) {
                accession = chr;
            }

            // HGVS genomic coordinate is 1-based; assume input position is already 1-based
            return accession + ":g." + position;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static boolean isValidHGVS(String notation) {
        boolean valid = notation != null && notation.matches("^[A-Z0-9_.]+:g\\.\\d+([ACGT]>[ACGT])?$");
        // IGV only supports genomic HGVS notation for now
        return valid && notation.contains(":g.");
    }

    /**
     * Searches for the given HGVS notation in the provided genome.
     * Returns a SearchResult with the corresponding chromosome and position if found,
     * otherwise returns null.
     */
    public static SearchCommand.SearchResult search(String hgvs, Genome genome) {
        if (!isValidHGVS(hgvs)) {
            return null;
        }

        String[] parts = hgvs.split(":g\\.");
        String accession = parts[0];
        String positionPart = parts[1];

        // Extract position (ignoring any variant information for now)
        String positionStr = positionPart.split("[ACGT]>[ACGT]")[0];
        int position = Integer.parseInt(positionStr);

        // Find chromosome corresponding to the accession
        String chr = genome.getCanonicalChrName(accession);

        // Return search result with chromosome and position
        return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, chr, position - 1, position);


    }

}