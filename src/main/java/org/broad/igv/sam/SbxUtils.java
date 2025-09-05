package org.broad.igv.sam;

import java.awt.*;
import java.util.HashSet;
import java.util.Set;

public class SbxUtils {
    public static final Color fullDiscordantInsertionColor = new Color(180, 180, 180);
    public static final Color discordantInsertionColor = new Color(140, 140, 140);
    public static final int DISCORDANT_BASE_QUALITY_MAX = 10; // Threshold for a base to be considered discordant

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
}
