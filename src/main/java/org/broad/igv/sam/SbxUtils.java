package org.broad.igv.sam;

import htsjdk.samtools.SAMRecord;
import org.broad.igv.sam.sbx.NullAlignment;

import java.awt.*;
import java.util.HashSet;
import java.util.Set;

public class SbxUtils {
    public static final Color fullDiscordantInsertionColor = new Color(180, 180, 180);
    public static final Color discordantInsertionColor = new Color(140, 140, 140);
    public static final int DISCORDANT_BASE_QUALITY_MAX = 10; // Threshold for a base to be considered discordant
    public static final int SBX_LOW_BASEQ_TAIL_THRESHOLD = 30;

    /*
     * Returns true if there is at least one discordant base within the insertion sequence
     */
    public static boolean isDiscordantInsertion(Alignment alignment, AlignmentBlock aBlock) {
        return aBlock != null && (isDiscordantInsertion(aBlock) || isDiscordantHpInsertion(alignment, aBlock) || hasAdjacentDiscordantBase(alignment, aBlock));
    }

    /*
     * Returns true if all bases in the insertion sequence are discordant
     */
    public static boolean isFullDiscordantInsertion(AlignmentBlock aBlock) {
        if (aBlock == null) {
            return false;
        }
        for (int i = 0; i < aBlock.getLength(); i++) {
            if (aBlock.getQuality(i) > DISCORDANT_BASE_QUALITY_MAX) {
                return false;
            }
        }
        return true;
    }

    /*
     * Returns true if there is at least one discordant base within the insertion sequence
     */
    public static boolean isDiscordantInsertion(AlignmentBlock aBlock) {
        for (int i = 0; i < aBlock.getLength(); i++) {
            if (aBlock.getQuality(i) <= DISCORDANT_BASE_QUALITY_MAX) {
                return true;
            }
        }
        return false;
    }

    /*
     * Returns true if the insertion is a homopolymer,and there is a discordant base either
     * within the insertion sequence or elsewhere in the same homopolymer sequence
     */
    public static boolean isDiscordantHpInsertion(Alignment alignment, AlignmentBlock aBlock) {
        byte base = aBlock.getBase(0);
        for (int i = 1; i < aBlock.getLength(); i++) {
            byte nextBase = aBlock.getBase(i);
            if (nextBase != base) {
                return false;
            }
        }
        int position = aBlock.getStart() - 1;
        while (position >= 0 && alignment.getBase(position) == base) {
            position--;
        }
        position++;
        while (position < alignment.getReadSequence().length() && alignment.getBase(position) == base) {
            if (alignment.getPhred(position) <= DISCORDANT_BASE_QUALITY_MAX) {
                return true;
            }
            position++;
        }
        return false;
    }

    /*
     * Returns true if one of the bases adjacent to the insertion block is discordant
     */
    public static boolean hasAdjacentDiscordantBase(Alignment alignment, AlignmentBlock aBlock) {
        int position = (int) aBlock.getStart(); // Genome position immediately after insertion
        if (alignment.getBase(position) != 0 && alignment.getPhred(position) <= DISCORDANT_BASE_QUALITY_MAX) {
            return true;
        }
        if (alignment.getBase(position - 1) != 0 && alignment.getPhred(position - 1) <= DISCORDANT_BASE_QUALITY_MAX) {
            return true;
        }
        return false;
    }

    /*
     * Generates a human-readable string of comma-delimited quality scores for displaying when clicking on an insertion
     */
    public static String qualityString(String rawQualString) {
        StringBuilder res = new StringBuilder();
        for (int i = 0; i < rawQualString.length(); i++) {
            if (i > 0) res.append(",");
            res.append((int) rawQualString.charAt(i));
        }
        return res.toString();
    }

    /**
     * Heuristic to determine if a set of quality scores is consistent with SBX alignments.
     *
     * @param qualityScores A set of quality scores observed in a collection of alignments
     * @return
     */
    public static boolean isSBXAlignments(Set<Integer> qualityScores) {

        Set<Integer> expectedScores1 = new HashSet<>(java.util.Arrays.asList(5, 22, 39));
        Set<Integer> expectedScores2 = new HashSet<>(java.util.Arrays.asList(0, 18, 93));

        if (qualityScores.size() == 0 || (qualityScores.size() == 1 && qualityScores.contains(0))) {
            return false;  // All zero quality scores, this can happen in non-SBX files
        }

        Set<Integer> expectedSet = null;
        for (Integer score : qualityScores) {

            if (expectedSet == null) {
                if (expectedScores1.contains(score)) {
                    expectedSet = expectedScores1;
                } else if (expectedScores2.contains(score)) {
                    expectedSet = expectedScores2;
                } else {
                    return false;  // Score not in either expected set
                }
            } else if (!expectedSet.contains(score)) {
                return false;  // Score not in the expected set
            }
        }

        return true;
    }

    public static int[] trimTails(Alignment alignment) {
        int simplexTailLeftEnd = 0;
        int alignmentLength = alignment.getEnd() - alignment.getStart();
        for (; simplexTailLeftEnd < alignmentLength && (alignment.getPhred(alignment.getStart() + simplexTailLeftEnd) <= SBX_LOW_BASEQ_TAIL_THRESHOLD || alignment.getBase(alignment.getStart() + simplexTailLeftEnd) == 0); simplexTailLeftEnd++) {

        }
        int simplexTailRightStart = alignmentLength - 1;
        for (; simplexTailRightStart >= 0 && (alignment.getPhred(alignment.getStart() + simplexTailRightStart) <= SBX_LOW_BASEQ_TAIL_THRESHOLD || alignment.getBase(alignment.getStart() + simplexTailRightStart) == 0); simplexTailRightStart--)
            ;

        return new int[]{simplexTailLeftEnd, simplexTailRightStart};
    }

    public static Alignment trimSimplexTails(SAMAlignment parent) {

        SAMRecord record = parent.getRecord().deepCopy();
        String cigarString = record.getCigarString();
        byte[] readBases = record.getReadBases();
        byte[] readBaseQualities = record.getBaseQualities();
        int start = 0;
        while (start < readBases.length && readBaseQualities[start] <= SBX_LOW_BASEQ_TAIL_THRESHOLD) {
            start++;
        }
        int end = readBases.length - 1;
        while (end >= 0 && readBaseQualities[end] <= SBX_LOW_BASEQ_TAIL_THRESHOLD) {
            end--;
        }

        if (start == 0 && end == readBases.length - 1) {
            // No trimming
            return parent;
        } else if (start > end) {
            // Entire read trimmed
            return NullAlignment.getInstance();
        }

        // Trim the read bases and qualities
        byte[] newReadBases = new byte[end - start + 1];
        byte[] newBaseQuals = new byte[end - start + 1];
        for (int i = 0; i < newReadBases.length; i++) {
            newReadBases[i] = readBases[start + i];
            newBaseQuals[i] = readBaseQualities[start + i];
        }
        record.setReadBases(newReadBases);
        record.setBaseQualities(newBaseQuals);

        // Trim the CIGAR string
        java.util.List<SAMAlignment.CigarOperator> operators = SAMAlignment.buildOperators(cigarString);
        StringBuilder newCigarString = new StringBuilder();

        int readStart = 0;
        int newAlignmentStart = record.getAlignmentStart();

        for (SAMAlignment.CigarOperator op : operators) {
            int length = op.nBases;
            char type = op.operator;
            final boolean consumesQuery = isConsumesQuery(type);
            final boolean consumesRef = isConsumesRef(type);

            int curEnd = readStart;
            if (consumesQuery) {
                curEnd += length;
            }

            if (curEnd <= start) {
                // Operation fully in left tail
                if (consumesRef) {
                    newAlignmentStart += length;
                }
            } else if (readStart > end) {
                // All done
                break;
            } else if (readStart >= start && curEnd <= end + 1) {
                // Operation fully retained
                newCigarString.append(length).append(type);
                readStart = curEnd;
            } else {
                // Operation partially retained
                int truncatedLeft = Math.max(0, start - readStart);
                int truncatedRight = Math.max(0, curEnd - (end + 1));
                int newLength = length - truncatedLeft - truncatedRight;
                if (newLength > 0) {
                    newCigarString.append(newLength).append(type);
                }
                if (consumesRef) {
                    newAlignmentStart += truncatedLeft;
                }
                readStart = curEnd;
            }
        }
        record.setAlignmentStart(newAlignmentStart);
        record.setCigarString(newCigarString.toString());
        return new SAMAlignment(record);
    }

    private static boolean isConsumesRef(char type) {
        boolean consumesRef = switch (type) {
            case SAMAlignment.DELETION, SAMAlignment.MISMATCH, SAMAlignment.MATCH,
                 SAMAlignment.PERFECT_MATCH, SAMAlignment.SKIPPED_REGION -> true;
            default -> false;
        };
        return consumesRef;
    }

    private static boolean isConsumesQuery(char type) {
        boolean consumesQuery = switch (type) {
            case SAMAlignment.INSERTION, SAMAlignment.MISMATCH, SAMAlignment.MATCH,
                 SAMAlignment.PERFECT_MATCH, SAMAlignment.SOFT_CLIP -> true;
            default -> false;
        };
        return consumesQuery;
    }

}