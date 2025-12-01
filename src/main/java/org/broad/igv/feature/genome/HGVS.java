// Java
package org.broad.igv.feature.genome;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.ui.action.SearchCommand;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class HGVS {

    public static boolean isValidHGVS(String notation) {
        if (notation == null) return false;
        // Accept genomic, coding, or protein HGVS with optional intronic offset and broad variant specifications.
        // Coding positions can be:
        //  - c.\d+      : CDS positions
        //  - c.-\d+     : 5' UTR positions upstream of the CDS start
        //  - c.*\d+     : 3' UTR positions downstream of the CDS end
        // Optional intronic offset (+N or -N) is allowed for c. forms.
        // Protein positions p.X (optionally with AA prefix/suffix, ranges, and simple variant suffixes)
        // Variants supported (genomic/coding/protein subsets):
        //  - Substitutions: N>N (letters beyond just ACGT)
        //  - del, delSEQ (e.g., delA, delAT)
        //  - insSEQ (e.g., insAT)
        //  - dup
        //  - delinsSEQ
        //  - frameshift (fs*) simplified acceptance
        //  - Optional ranges X_Y before the variant
        String codingBase = "(?:\\d+|-\\d+|\\*\\d+)";
        String codingRange = codingBase + "(?:_" + codingBase + ")?"; // allow optional range X_Y
        String intronOffset = "(?:[+-]\\d+)?";
        String subst = "[A-Za-z]+>[A-Za-z]+";
        String del = "del(?:[A-Za-z]+)?"; // allow delA, delAT, or just del
        String ins = "ins[A-Za-z]+"; // require sequence for insertions
        String dup = "dup";
        String delins = "delins[A-Za-z]+";
        String fs = "fs\\*?\\d*"; // frameshift notation
        String variant = "(?:" + subst + "|" + delins + "|" + ins + "|" + del + "|" + dup + "|" + fs + ")?";
        String genomic = "g\\.\\d+" + "(?:" + subst + "|" + delins + "|" + ins + "|" + del + "|" + dup + ")?"; // genomic doesn't use fs
        String coding = "c\\." + codingRange + intronOffset + variant;
        // Protein: optional AA prefix (1-3 letters or *), position, optional range, optional AA suffix or variant token
        String aa = "[A-Za-z*]{0,3}";
        String proteinPos = aa + "\\d+" + "(?:_" + aa + "\\d+)?"; // allow range with AA context
        // Protein variants can have AA names before fs (e.g., Thr758Serfs or Leu29Alafs*18)
        String proteinVariant = "(?:" + aa + "(?:" + fs + "|" + del + "|" + ins + "|" + dup + "|" + delins + ")?|" + del + "|" + ins + "|" + dup + "|" + delins + "|" + fs + ")?";
        String protein = "p\\." + proteinPos + proteinVariant;
        return notation.matches("^[A-Za-z0-9_.]+:(?:" + genomic + "|" + coding + "|" + protein + ")$");
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

        // Determine type and extract accession and position
        int idxG = hgvs.indexOf(":g.");
        int idxC = hgvs.indexOf(":c.");
        int idxP = hgvs.indexOf(":p.");
        String type;
        int idx;
        if (idxG >= 0) {
            type = "g";
            idx = idxG;
        } else if (idxC >= 0) {
            type = "c";
            idx = idxC;
        } else if (idxP >= 0) {
            type = "p";
            idx = idxP;
        } else {
            return null;
        }
        String accession = hgvs.substring(0, idx);
        String positionPart = hgvs.substring(idx + 3); // skip ':g.' or ':c.' or ':p.'

        if ("g".equals(type)) {
            if (positionPart == null) return null;
            Matcher matcher = Pattern.compile("^(\\d+)").matcher(positionPart);
            if (!matcher.find()) return null;
            int position = Integer.parseInt(matcher.group(1));
            String chr = genome.getCanonicalChrName(accession);
            return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, chr, position - 1, position);
        } else if ("p".equals(type)) {
            // Protein position mapping: map codon(s) to genomic span.
            final BasicFeature transcript = getTranscript(genome, accession);
            if (transcript == null) return null;

            // Strip trailing variant info after first non position pattern
            String proteinPart = positionPart; // e.g. A123T, 123_125del, Gly123_Val125delins...
            // Regex to capture first and optional second numeric protein positions
            Matcher pm = Pattern.compile("^[A-Za-z*]{0,3}(\\d+)(?:_[A-Za-z*]{0,3}(\\d+))?").matcher(proteinPart);
            if (!pm.find()) return null;
            int p1 = Integer.parseInt(pm.group(1));
            String p2Str = pm.group(2);
            int p2 = p1;
            if (p2Str != null) {
                p2 = Integer.parseInt(p2Str);
                // Strand: negative strand higher protein position => lower genomic coordinate, choose strand-aware start/end later
            }
            // Get codon(s)
            var codon1 = transcript.getCodon(genome, transcript.getChr(), p1);
            if (codon1 == null || !codon1.isGenomePositionsSet()) return null;
            int start1 = Integer.MAX_VALUE;
            int end1 = Integer.MIN_VALUE;
            for (int gp : codon1.getGenomePositions()) {
                start1 = Math.min(start1, gp);
                end1 = Math.max(end1, gp);
            }
            int regionStart = start1;
            int regionEnd = end1; // inclusive base
            if (p2 != p1) {
                var codon2 = transcript.getCodon(genome, transcript.getChr(), p2);
                if (codon2 == null || !codon2.isGenomePositionsSet()) return null;
                int start2 = Integer.MAX_VALUE;
                int end2 = Integer.MIN_VALUE;
                for (int gp : codon2.getGenomePositions()) {
                    start2 = Math.min(start2, gp);
                    end2 = Math.max(end2, gp);
                }
                if (transcript.getStrand() == Strand.POSITIVE) {
                    regionStart = Math.min(start1, start2);
                    regionEnd = Math.max(end1, end2);
                } else {
                    // Negative strand: lower genomic coordinate is further downstream; choose min for start
                    regionStart = Math.min(start1, start2);
                    regionEnd = Math.max(end1, end2);
                }
            }
            // Convert inclusive end to half-open end for SearchResult (end position is 1-based exclusive)
            int halfOpenEnd = regionEnd + 1;
            return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), regionStart, halfOpenEnd);
        } else {

            // Coding / UTR mapping
            BasicFeature transcript = getTranscript(genome, accession);
            if (transcript != null) {
                // UTR 5' c.-N with optional range and intronic offset (e.g., c.-211_-215 or c.-211-1058C>G)
                Matcher utr5Matcher = Pattern.compile("^-(\\d+)(?:_-(\\d+))?([+-]\\d+)?").matcher(positionPart);
                if (utr5Matcher.find()) {
                    int n = Integer.parseInt(utr5Matcher.group(1));
                    // For ranges, use the second position if it exists
                    String n2Str = utr5Matcher.group(2);
                    if (n2Str != null) {
                        int n2 = Integer.parseInt(n2Str);
                        // Use the higher c.- number to get the genomic start position
                        n = n2;
                    }
                    int firstCodingGenomic = transcript.codingToGenomePosition(1);
                    if (firstCodingGenomic > 0) {
                        int pos = transcript.getStrand() == Strand.POSITIVE ? (firstCodingGenomic - n) : (firstCodingGenomic + n);
                        // Apply intronic offset if present
                        String offsetStr = utr5Matcher.group(3);
                        if (offsetStr != null) {
                            int offset = Integer.parseInt(offsetStr);
                            if (transcript.getStrand() == Strand.NEGATIVE) offset = -offset;
                            pos += offset;
                        }
                        return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), pos - 1, pos);
                    }
                    return null;
                }

                // UTR 3' c.*N with optional range and intronic offset (e.g., c.*526_*529delATCA or c.*123+45)
                Matcher utr3Matcher = Pattern.compile("^\\*(\\d+)(?:_\\*(\\d+))?([+-]\\d+)?").matcher(positionPart);
                if (utr3Matcher.find()) {
                    int n = Integer.parseInt(utr3Matcher.group(1));
                    // For ranges, use the second position if it exists
                    // On negative strand, higher c.* numbers map to lower genomic coordinates,
                    // so we want the larger number for the genomic start
                    String n2Str = utr3Matcher.group(2);
                    if (n2Str != null) {
                        int n2 = Integer.parseInt(n2Str);
                        // For negative strand, use the higher c.* position (which maps to lower genomic coord)
                        // For positive strand, use the lower c.* position (which maps to lower genomic coord)
                        if (transcript.getStrand() == Strand.NEGATIVE) {
                            n = n2; // Use the second (higher) position for negative strand
                        }
                        // For positive strand, keep n as is (the first, lower position)
                    }
                    int codingLen = 0;
                    if (transcript.getExons() != null) {
                        for (var exon : transcript.getExons()) {
                            codingLen += exon.getCodingLength();
                        }
                    }
                    if (codingLen > 0) {
                        int lastCodingGenomic = transcript.codingToGenomePosition(codingLen);
                        if (lastCodingGenomic > 0) {
                            int pos = transcript.getStrand() == Strand.POSITIVE ? (lastCodingGenomic + n) : (lastCodingGenomic - n);
                            // Apply intronic offset if present
                            String offsetStr = utr3Matcher.group(3);
                            if (offsetStr != null) {
                                int offset = Integer.parseInt(offsetStr);
                                if (transcript.getStrand() == Strand.NEGATIVE) offset = -offset;
                                pos += offset;
                            }
                            return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), pos - 1, pos);
                        }
                    }
                    return null;
                }

                // CDS position with optional range and intronic offset (e.g., c.123_124insAT or c.505-2A>G)
                Matcher matcher = Pattern.compile("^(\\d+)(?:_(\\d+))?([+-]\\d+)?").matcher(positionPart);
                if (!matcher.find()) return null;
                int codingPosition = Integer.parseInt(matcher.group(1));
                // For ranges in CDS, use the first position for positive strand, second for negative
                String pos2Str = matcher.group(2);
                if (pos2Str != null) {
                    int codingPosition2 = Integer.parseInt(pos2Str);
                    // On negative strand, higher coding positions map to lower genomic coordinates
                    if (transcript.getStrand() == Strand.NEGATIVE) {
                        codingPosition = codingPosition2; // Use the higher coding position
                    }
                    // On positive strand, keep the lower coding position
                }
                int genomicPosition = transcript.codingToGenomePosition(codingPosition);
                String offsetStr = matcher.group(3);
                if (offsetStr != null) {
                    int offset = Integer.parseInt(offsetStr);
                    if (transcript.getStrand() == Strand.NEGATIVE) offset = -offset;
                    genomicPosition += offset;
                }
                if (genomicPosition > 0) {
                    return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), genomicPosition - 1, genomicPosition);
                }
            }
            return null;
        }
    }

    private static BasicFeature getTranscript(Genome genome, String accession) {

        IGVFeature feature = genome.getManeFeature(accession);

        if (feature == null && genome.getFeatureDB() != null) {
            var nf = genome.getFeatureDB().getFeature(accession);
            if (nf instanceof IGVFeature) {
                feature = (IGVFeature) nf;
            }
        }

        if (!(feature instanceof BasicFeature)) {
            return null;
        }
        BasicFeature transcript = (BasicFeature) feature;
        return transcript;
    }

    /**
     * Returns genomic HGVS notation: <RefSeqAccession>:g.<position>
     * Example: NC_000001.11:g.1234567
     */
    public static String getHGVSPosition(Genome genome, String chr, int position) {
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


}