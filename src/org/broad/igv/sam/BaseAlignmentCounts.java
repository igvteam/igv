/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date 11/29/11
 */
abstract public class BaseAlignmentCounts implements AlignmentCounts {

    private static Logger log = Logger.getLogger(BaseAlignmentCounts.class);

    private static char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    private static Map<String, Set<Integer>> knownSnps;
    int start;
    int end;
    protected boolean countDeletedBasesCovered = false;


    private BisulfiteCounts bisulfiteCounts;


    public BaseAlignmentCounts(int start, int end, AlignmentTrack.BisulfiteContext bisulfiteContext) {
        final PreferenceManager prefs = PreferenceManager.getInstance();
        String snpsFile = prefs.get(PreferenceManager.KNOWN_SNPS, null);
        if (snpsFile != null && knownSnps == null) {
            loadKnownSnps(snpsFile);
        }
        this.start = start;
        this.end = end;

        countDeletedBasesCovered = prefs.getAsBoolean(PreferenceManager.SAM_COUNT_DELETED_BASES_COVERED);

        if (!Globals.isHeadless() && bisulfiteContext != null) {
            bisulfiteCounts = new BisulfiteCounts(bisulfiteContext, GenomeManager.getInstance().getCurrentGenome());
        }

    }


    public int getStart() {
        return start;
    }


    public int getEnd() {
        return end;
    }

    public String getChr() {
        return null;
    }


    public BisulfiteCounts getBisulfiteCounts() {
        return bisulfiteCounts;
    }

    /**
     * Increment the counts for this alignment.   Does not consider softclips.
     *
     * @param alignment
     */
    public void incCounts(Alignment alignment) {

        if (bisulfiteCounts != null) {
            bisulfiteCounts.incrementCounts(alignment);
        }

        int alignmentStart = alignment.getAlignmentStart();
        int alignmentEnd = alignment.getAlignmentEnd();

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        if (blocks != null) {
            int lastBlockEnd = -1;
            int gapIdx = 0;
            char[] gapTypes = alignment.getGapTypes();

            for (AlignmentBlock b : blocks) {
                if (b.getEnd() < start) continue;
                if (b.getStart() > end) break;

                //Strand strand = alignment.getFirstOfPairStrand();
                Strand strand = alignment.getReadStrand();

                // Don't count softclips
                if (!b.isSoftClipped() && strand != Strand.NONE) {
                    final boolean isNegativeStrand = strand == Strand.NEGATIVE;

                    incBlockCounts(b, isNegativeStrand); //alignment.isNegativeStrand());

                    // Count deletions
                    if (gapTypes != null
                            && (lastBlockEnd >= 0)
                            && (gapIdx < gapTypes.length)) {
                        if (gapTypes[gapIdx] == SamAlignment.DELETION) {
                            for (int pos = lastBlockEnd; pos < b.getStart(); pos++) {
                                incrementDeletion(pos, isNegativeStrand);
                            }
                        }
                        gapIdx++;
                    }
                    lastBlockEnd = b.getEnd();
                }
            }
            // Count insertions
            final AlignmentBlock[] insertions = alignment.getInsertions();
            if (insertions != null) {
                for (AlignmentBlock insBlock : insertions) {
                    if (insBlock.getEnd() < start) continue;
                    if (insBlock.getStart() > end) break;

                    incrementInsertion(insBlock);
                }
            }
        } else {
            // No alignment blocks => no bases (e.g. .aln or .aligned files).  Just do total count.
            for (int pos = alignmentStart; pos < alignmentEnd; pos++) {
                byte q = 0;
                incPositionCount(pos, (byte) 'n', q, alignment.isNegativeStrand());
            }
        }
    }

    public String getValueStringAt(int pos) {

        StringBuffer buf = new StringBuffer();
        int totalCount = getTotalCount(pos);
        buf.append("Total count: " + totalCount);
        for (char c : nucleotides) {
            int negCount = getNegCount(pos, (byte) c);
            int posCount = getPosCount(pos, (byte) c);
            int count = negCount + posCount;
            int percent = (int) Math.round(((float) count) * 100 / totalCount);
            char cU = Character.toUpperCase(c);
            buf.append("<br>" + cU + "      : " + count);
            if (count != 0) {
                buf.append("  (" + percent + "%,     " + posCount + "+,   " + negCount + "- )");
            }
        }

        int delCount = getDelCount(pos);
        if (delCount > 0) {
            buf.append("<br>DEL: " + delCount);
        }
        int insCount = getInsCount(pos);
        if (insCount > 0) {
            buf.append("<br>INS: " + insCount);
        }

        return buf.toString();

    }

    public boolean isMismatch(int pos, byte ref, String chr, float snpThreshold) {

        Set<Integer> filteredSnps = knownSnps == null ? null : knownSnps.get(chr);
        if (filteredSnps == null || !filteredSnps.contains(pos + 1)) {
            float threshold = snpThreshold * getTotalQuality(pos);
            if (ref > 0) {
                if (ref < 96) ref += 32;  // a fast "toLowercase"
                for (char c : nucleotides) {
                    if (c != ref && c != 'n' && getQuality(pos, (byte) c) > threshold) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Load the set of known snps from a tab delimited file, format
     * chr < tab> location
     * The location is "1 base"  (first nucleotide is position 1).
     *
     * @param snpFile
     */
    private static synchronized void loadKnownSnps(String snpFile) {

        // This method might get called many times concurrently, but we only want to load these once.
        if (knownSnps != null) {
            return;
        }

        knownSnps = new HashMap();
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(new ResourceLocator(snpFile));
            String nextLine = "";
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0];
                Set<Integer> snps = knownSnps.get(chr);
                if (snps == null) {
                    snps = new HashSet(10000);
                    knownSnps.put(chr, snps);
                }
                snps.add(new Integer(tokens[1]));
            }
        } catch (Exception e) {
            knownSnps = null;
            log.error("", e);
            MessageUtils.showMessage("Error loading snps file: " + snpFile + " (" + e.toString() + ")");
        } finally {
            reader.close();
        }


    }

    protected abstract void incPositionCount(int pos, byte n, byte q, boolean negativeStrand);

    protected abstract void incrementInsertion(AlignmentBlock insBlock);

    protected abstract void incrementDeletion(int pos, boolean negativeStrand);

    protected abstract void incBlockCounts(AlignmentBlock b, boolean b1);


}
