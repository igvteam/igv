// Java
package org.igv.feature.genome;

import htsjdk.samtools.util.SequenceUtil;
import org.igv.feature.BasicFeature;
import org.igv.feature.Exon;
import org.igv.feature.IGVFeature;
import org.igv.feature.Strand;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.action.SearchCommand;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class HGVS {

    private static final Logger log = LogManager.getLogger(HGVS.class);

    public static boolean isValidHGVS(String notation) {
        if (notation == null) return false;
        // We only need to validate that we can parse the notation in the search method.
        // Check for basic structure: <accession>(<gene>)?:g.|c.|n.|p. followed by a position we can parse.
        // We don't validate the variant details.

        // Genomic: g.\d+ (with optional range and anything after)
        String genomic = "g\\.\\d+.*";
        // Coding: c. followed by optional -, *, then digits, with optional intronic offset and anything after
        String coding = "c\\.[-*]?\\d+.*";
        // Non-coding: n. followed by optional leading '-' then digits, anything after
        String nonCoding = "n\\.-?\\d+.*";
        // Protein: p. followed by optional AA letters, digits, with optional range and anything after
        String protein = "p\\.[A-Za-z*]*\\d+.*";

        // Optional gene symbol in parentheses immediately after accession
        String accessionWithOptionalGene = "^[A-Za-z0-9_.]+(?:\\([^)]+\\))?";

        return notation.matches(accessionWithOptionalGene + ":(?:" + genomic + "|" + coding + "|" + nonCoding + "|" + protein + ")$");
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
        int idxN = hgvs.indexOf(":n.");
        int idxP = hgvs.indexOf(":p.");
        String type;
        int idx;
        if (idxG >= 0) {
            type = "g";
            idx = idxG;
        } else if (idxC >= 0) {
            type = "c";
            idx = idxC;
        } else if (idxN >= 0) {
            type = "n";
            idx = idxN;
        } else if (idxP >= 0) {
            type = "p";
            idx = idxP;
        } else {
            return null;
        }
        String accession = hgvs.substring(0, idx);
        // Strip optional trailing gene symbol in parentheses, e.g., "NM_000302.3(PLOD1)" -> "NM_000302.3"
        if (accession.endsWith(")")) {
            int openIdx = accession.lastIndexOf('(');
            if (openIdx > 0) {
                accession = accession.substring(0, openIdx);
            }
        }
        String positionPart = hgvs.substring(idx + 3); // skip ':g.' or ':c.' or ':n.' or ':p.'

        switch (type) {
            case "g": {
                // Handle both single positions (123) and ranges (123_456)
                Matcher matcher = Pattern.compile("^(\\d+)(?:_(\\d+))?").matcher(positionPart);
                if (!matcher.find()) return null;
                int start = Integer.parseInt(matcher.group(1));
                String endGroup = matcher.group(2);
                int end = endGroup != null ? Integer.parseInt(endGroup) : start;
                String chr = genome.getCanonicalChrName(accession);
                // UCSC style: 0-based start, half-open end. HGVS g. is 1-based inclusive.
                return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, chr, start - 1, end);
            }
            case "p": {
                // Protein position mapping: map codon(s) to genomic span.
                final BasicFeature transcript = getTranscript(genome, accession);
                if (transcript == null) return null;

                // Regex to capture first and optional second numeric protein positions (skip AA labels)
                Matcher pm = Pattern.compile("^[A-Za-z*]{0,3}(\\d+)(?:_[A-Za-z*]{0,3}(\\d+))?").matcher(positionPart);
                if (!pm.find()) return null;
                int p1 = Integer.parseInt(pm.group(1));
                String p2Str = pm.group(2);
                int p2 = p2Str != null ? Integer.parseInt(p2Str) : p1;
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
                    regionStart = Math.min(start1, start2);
                    regionEnd = Math.max(end1, end2);
                }
                // Convert inclusive end to half-open end for SearchResult (end position is 1-based exclusive)
                int halfOpenEnd = regionEnd + 1;
                return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), regionStart, halfOpenEnd);
            }
            case "n": {
                // Non-coding transcript mapping: n.123 or n.-123 maps relative to transcript start
                BasicFeature transcript = getTranscript(genome, accession);
                if (transcript == null) return null;

                // Parse signed position with optional range and intronic offset (e.g., n.123, n.123_456, n.-7080_-1781, n.123+5)
                Matcher matcher = Pattern.compile("^(-?\\d+)(?:_(-?\\d+))?([+-]\\d+)?").matcher(positionPart);
                if (!matcher.find()) return null;

                int t1 = Integer.parseInt(matcher.group(1));
                String t2Str = matcher.group(2);

                // Determine transcript positions for both ends; if single, both are the same
                int t2 = (t2Str != null) ? Integer.parseInt(t2Str) : t1;

                // Map both transcript positions to genomic
                int g1 = transcriptPositionToGenomicPosition(transcript, t1);
                int g2 = transcriptPositionToGenomicPosition(transcript, t2);
                if (g1 <= 0 || g2 <= 0) return null;

                // Apply intronic offset (if any) to BOTH endpoints, strand-aware
                String offsetStr = matcher.group(3);
                if (offsetStr != null) {
                    int offset = Integer.parseInt(offsetStr);
                    if (transcript.getStrand() == Strand.NEGATIVE) offset = -offset;
                    g1 += offset;
                    g2 += offset;
                }

                // Normalize to genomic span regardless of strand
                int regionStart = Math.min(g1, g2);
                int regionEndInclusive = Math.max(g1, g2);
                int halfOpenEnd = regionEndInclusive + 1;
                return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), regionStart, halfOpenEnd);
            }
            case "c": {
                // Coding / UTR mapping
                BasicFeature transcript = getTranscript(genome, accession);
                if (transcript != null) {
                    // UTR 5' c.-N with optional range and intronic offset (e.g., c.-211_-215 or c.-211-1058C>G)
                    Matcher utr5Matcher = Pattern.compile("^-(\\d+)(?:_-(\\d+))?([+-]\\d+)?").matcher(positionPart);
                    if (utr5Matcher.find()) {
                        int n1 = Integer.parseInt(utr5Matcher.group(1));
                        String n2Str = utr5Matcher.group(2);
                        Integer n2 = n2Str != null ? Integer.parseInt(n2Str) : null;
                        int firstCodingGenomic = transcript.codingToGenomePosition(1);
                        if (firstCodingGenomic > 0) {
                            int g1 = transcript.getStrand() == Strand.POSITIVE ? (firstCodingGenomic - n1) : (firstCodingGenomic + n1);
                            int g2 = g1;
                            if (n2 != null) {
                                g2 = transcript.getStrand() == Strand.POSITIVE ? (firstCodingGenomic - n2) : (firstCodingGenomic + n2);
                            }
                            // Apply intronic offset (single value) to both ends if present
                            String offsetStr = utr5Matcher.group(3);
                            if (offsetStr != null) {
                                int offset = Integer.parseInt(offsetStr);
                                if (transcript.getStrand() == Strand.NEGATIVE) offset = -offset;
                                g1 += offset;
                                g2 += offset;
                            }
                            int start = Math.min(g1, g2);
                            int endInclusive = Math.max(g1, g2);
                            int endExclusive = endInclusive + 1;
                            return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), start, endExclusive);
                        }
                        return null;
                    }

                    // UTR 3' c.*N with optional range and intronic offset (e.g., c.*526_*529delATCA or c.*123+45)
                    Matcher utr3Matcher = Pattern.compile("^\\*(\\d+)(?:_\\*(\\d+))?([+-]\\d+)?").matcher(positionPart);
                    if (utr3Matcher.find()) {
                        int n1 = Integer.parseInt(utr3Matcher.group(1));
                        String n2Str = utr3Matcher.group(2);
                        Integer n2 = n2Str != null ? Integer.parseInt(n2Str) : null;
                        int codingLen = 0;
                        if (transcript.getExons() != null) {
                            for (var exon : transcript.getExons()) {
                                codingLen += exon.getCodingLength();
                            }
                        }
                        if (codingLen > 0) {
                            int lastCodingGenomic = transcript.codingToGenomePosition(codingLen);
                            if (lastCodingGenomic > 0) {
                                int g1 = transcript.getStrand() == Strand.POSITIVE ? (lastCodingGenomic + n1) : (lastCodingGenomic - n1);
                                int g2 = g1;
                                if (n2 != null) {
                                    g2 = transcript.getStrand() == Strand.POSITIVE ? (lastCodingGenomic + n2) : (lastCodingGenomic - n2);
                                }
                                // Apply intronic offset (single value) to both ends if present
                                String offsetStr = utr3Matcher.group(3);
                                if (offsetStr != null) {
                                    int offset = Integer.parseInt(offsetStr);
                                    if (transcript.getStrand() == Strand.NEGATIVE) offset = -offset;
                                    g1 += offset;
                                    g2 += offset;
                                }
                                int start = Math.min(g1, g2);
                                int endInclusive = Math.max(g1, g2);
                                int endExclusive = endInclusive + 1;
                                return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), start, endExclusive);
                            }
                        }
                        return null;
                    }

                    // CDS position with optional range
                    // First parse endpoints c.X(_Y)? ignoring intronic offsets
                    Matcher cpos = Pattern.compile("^(\\d+)(?:_(\\d+))?").matcher(positionPart);
                    if (!cpos.find()) return null;
                    int c1 = Integer.parseInt(cpos.group(1));
                    String c2Str = cpos.group(2);
                    int c2 = c2Str != null ? Integer.parseInt(c2Str) : c1;

                    // Map both coding positions to genomic
                    int g1 = transcript.codingToGenomePosition(c1);
                    int g2 = transcript.codingToGenomePosition(c2);
                    if (g1 <= 0 || g2 <= 0) return null;

                    // Now parse optional intronic offsets for each endpoint separately
                    // Patterns like: 123+5 or 123-2 at the beginning, optionally followed by _ and second with offset
                    Matcher offs = Pattern.compile("^(\\d+)([+-]\\d+)?(?:_(\\d+)([+-]\\d+)?)?").matcher(positionPart);
                    if (offs.find()) {
                        String off1Str = offs.group(2);
                        String off2Str = offs.group(4);
                        if (off1Str != null) {
                            int off1 = Integer.parseInt(off1Str);
                            if (transcript.getStrand() == Strand.NEGATIVE) off1 = -off1;
                            g1 += off1;
                        }
                        if (off2Str != null) {
                            int off2 = Integer.parseInt(off2Str);
                            if (transcript.getStrand() == Strand.NEGATIVE) off2 = -off2;
                            g2 += off2;
                        }
                    }

                    // If there is no explicit second coding position, ensure single-site locus
                    if (c2Str == null) {
                        g2 = g1;
                    }

                    int start = Math.min(g1, g2);
                    int endInclusive = Math.max(g1, g2);
                    int endExclusive = endInclusive + 1;
                    return new SearchCommand.SearchResult(SearchCommand.ResultType.LOCUS, transcript.getChr(), start, endExclusive);
                }
                return null;
            }
            default:
                return null;
        }
    }

    /**
     * Convert a transcript position (1-based, from transcription start) to genomic position
     * for non-coding transcripts. Walks through exons to find the genomic coordinate.
     */
    private static int transcriptPositionToGenomicPosition(BasicFeature transcript, int transcriptPos) {
        // Handle positions upstream of transcript start (negative n. values)
        if (transcriptPos <= 0) {
            int d = Math.abs(transcriptPos);
            return transcript.getStrand() == Strand.POSITIVE ? (transcript.getStart() - d) : (transcript.getEnd() + d);
        }

        List<Exon> exons = transcript.getExons();
        if (exons == null || exons.isEmpty()) {
            // No exons, treat as simple feature
            if (transcript.getStrand() == Strand.POSITIVE) {
                return transcript.getStart() + transcriptPos - 1;
            } else {
                return transcript.getEnd() - transcriptPos + 1;
            }
        }

        boolean positive = transcript.getStrand() == Strand.POSITIVE;
        int accumulatedLength = 0;

        // Sort exons appropriately based on strand
        List<Exon> sortedExons = new ArrayList<>(exons);
        if (!positive) {
            sortedExons.sort((e1, e2) -> Integer.compare(e2.getStart(), e1.getStart()));
        } else {
            sortedExons.sort(Comparator.comparingInt(Exon::getStart));
        }

        for (Exon exon : sortedExons) {
            int exonLength = exon.getEnd() - exon.getStart();
            if (accumulatedLength + exonLength >= transcriptPos) {
                // Position is in this exon
                int offsetInExon = transcriptPos - accumulatedLength - 1;
                if (positive) {
                    return exon.getStart() + offsetInExon;
                } else {
                    return exon.getEnd() - offsetInExon - 1;
                }
            }
            accumulatedLength += exonLength;
        }

        // Position beyond transcript end
        return -1;
    }

    private static BasicFeature getTranscript(Genome genome, String accession) {

        IGVFeature feature = genome.getManeTranscript(accession);

        if (feature == null && genome.getFeatureDB() != null) {
            var nf = genome.getFeatureDB().getFeature(accession);
            if (nf instanceof IGVFeature) {
                feature = (IGVFeature) nf;
            }
        }

        if (!(feature instanceof BasicFeature)) {
            return null;
        }
        return (BasicFeature) feature;
    }

    /**
     * Returns HGVS annotation for the position, for ref and alt bases.  If a MANE transcript is available that is
     * used with coding notation (c.), otherwise genome position is used with genomic notation (g.).
     * Example: NM_000302.3:c.1234 or NM_000302.3:c.123+5 (intronic) or NC_000001.11:g.1234567
     *
     * @param genome The genome
     * @param chr The chromosome name
     * @param position The genomic position (0-based)
     * @return HGVS notation string, or null if error
     */
    public static String createHGVSAnnotation(Genome genome, String chr, int position, byte reference, byte alternate) {

        try {
            BasicFeature transcript = genome.getManeTranscriptAt(chr, position);

            if (transcript != null && transcript.getExons() != null) {

                if (transcript.getStrand() == Strand.NEGATIVE) {
                    reference = SequenceUtil.complement(reference);
                    alternate = SequenceUtil.complement(alternate);
                }

                String positionString = "";

                String transcriptName = transcript.getName();
                for (String a : transcript.getAttributes().values()) {
                    if (a.startsWith("NM_") || a.startsWith("NR_")) {
                        transcriptName = a;
                        break;
                    }
                }
                if (transcriptName != null && !transcriptName.isEmpty()) {
                    // Check if position is within an exon (coding or non-coding)
                    boolean positionIsInExon = false;
                    for (Exon exon : transcript.getExons()) {
                        if (position >= exon.getStart() && position < exon.getEnd()) {
                            positionIsInExon = true;
                            break;
                        }
                    }

                    if (positionIsInExon) {
                        // Try to convert to coding position
                        int codingPosition = transcript.genomeToCodingPosition(position);

                        if (codingPosition >= 0) {
                            // Position is in a coding region, return c. notation (1-based)
                            positionString = transcriptName + ":c." + (codingPosition + 1);
                        } else {
                            // Position is in an exon but not coding - check if in UTR
                            int firstCodingPos = transcript.codingToGenomePosition(1);
                            if (firstCodingPos > 0) {
                                // Calculate total coding length
                                int codingLen = 0;
                                for (Exon exon : transcript.getExons()) {
                                    codingLen += exon.getCodingLength();
                                }
                                int lastCodingPos = transcript.codingToGenomePosition(codingLen);

                                boolean positive = transcript.getStrand() == Strand.POSITIVE;

                                // Check if in 5' UTR
                                if ((positive && position < firstCodingPos) || (!positive && position > firstCodingPos)) {
                                    int distance = Math.abs(position - firstCodingPos);
                                    positionString = transcriptName + ":c.-" + distance;
                                }
                                // Check if in 3' UTR
                                else if ((positive && position >= lastCodingPos) || (!positive && position <= lastCodingPos)) {
                                    int distance = Math.abs(position - lastCodingPos) + 1;
                                    positionString = transcriptName + ":c.*" + distance;
                                }
                            }
                        }
                    } else {
                        // Position is intronic - find nearest exon boundary
                        // For HGVS, we reference the last coding base in the nearest exon
                        int nearestExonEdge = -1;
                        int nearestCodingPos = -1;
                        int minDistance = Integer.MAX_VALUE;
                        boolean positive = transcript.getStrand() == Strand.POSITIVE;

                        for (Exon exon : transcript.getExons()) {
                            if (exon.getCodingLength() == 0) continue; // Skip non-coding exons

                            // Check distance to the last coding base at the start side of the exon
                            // exon.getStart() is 0-based inclusive
                            int distToStart = Math.abs(position - exon.getStart());
                            if (distToStart > 0 && distToStart < minDistance) {
                                minDistance = distToStart;
                                nearestExonEdge = exon.getStart();
                                // Get coding position of first base in this exon
                                nearestCodingPos = transcript.genomeToCodingPosition(exon.getCdStart());
                            }

                            // Check distance to the last coding base at the end side of the exon
                            // exon.getEnd() is 0-based exclusive, so last base is at getEnd()-1
                            int distToEnd = Math.abs(position - (exon.getEnd() - 1));
                            if (distToEnd > 0 && distToEnd < minDistance) {
                                minDistance = distToEnd;
                                nearestExonEdge = exon.getEnd() - 1;
                                // Get coding position of last base in this exon
                                nearestCodingPos = transcript.genomeToCodingPosition(exon.getCdEnd() - 1);
                            }
                        }

                        if (nearestCodingPos >= 0) {
                            // Calculate offset: positive = downstream of exon, negative = upstream of exon
                            int offset = position - nearestExonEdge;
                            // For positive strand: + means to the right, - means to the left
                            // For negative strand: + means to the left (genomically), - means to the right
                            // But in HGVS, the sign is relative to transcript direction, so we need to flip for negative strand
                            if (!positive) {
                                offset = -offset;
                            }
                            String sign = offset >= 0 ? "+" : "";
                            positionString = transcriptName + ":c." + (nearestCodingPos + 1) + sign + offset;
                        }
                    }
                }

                return positionString + (char) reference + ">" + (char) alternate;
            }

            // Fallback to genomic notation
            ChromAlias aliasRecord = genome.getAliasRecord(chr);
            String accession = chr;

            if (aliasRecord != null) {
                for (String alias : aliasRecord.values()) {
                    if (alias.startsWith("NC_") || alias.startsWith("NT_") || alias.startsWith("NW_") ||
                            alias.startsWith("NG_") || alias.startsWith("NM_") || alias.startsWith("NR_") ||
                            alias.startsWith("NP_")) {
                        accession = alias;
                        break;
                    }
                }
            }


            // HGVS genomic coordinate is 1-based; position parameter is 0-based
            return accession + ":g." + (position + 1) + (char) reference + ">" + (char) alternate;
        } catch (IOException e) {
            log.error("Error getting HGVS position", e);
            return null;
        }
    }


}
