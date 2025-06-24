/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.sam;

import htsjdk.tribble.readers.AsciiLineReader;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.mods.BaseModificationCounts;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author Jim Robinson
 * @date 11/29/11
 */
abstract public class BaseAlignmentCounts implements AlignmentCounts {

    private static Logger log = LogManager.getLogger(BaseAlignmentCounts.class);

    public static final char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    private static Map<String, Set<Integer>> knownSnps;
    int start;
    int end;
    private int bucketSize = 1;
    protected boolean countDeletedBasesCovered = false;

    private BisulfiteCounts bisulfiteCounts;
    private BaseModificationCounts baseModificationCounts;


    public BaseAlignmentCounts(int start, int end, AlignmentTrack.BisulfiteContext bisulfiteContext) {
        final IGVPreferences prefs = PreferencesManager.getPreferences();
        String snpsFile = prefs.get(KNOWN_SNPS, null);
        if (snpsFile != null && knownSnps == null) {
            loadKnownSnps(snpsFile);
        }
        this.start = start;
        this.end = end;

        countDeletedBasesCovered = prefs.getAsBoolean(SAM_COUNT_DELETED_BASES_COVERED);

        if (!Globals.isHeadless() && bisulfiteContext != null) {
            bisulfiteCounts = new BisulfiteCounts(bisulfiteContext, GenomeManager.getInstance().getCurrentGenome());
        }

        baseModificationCounts = new BaseModificationCounts();

    }

    /**
     * Get the set of base characters that are observed in this pileup.  This might characters other than 'a', 'c',
     * 'g', 't', and 'n', such as ambiguity codes or '='.
     * @return
     */
    public abstract Set<Byte> getBases();

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String getChr() {
        return null;
    }

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public BisulfiteCounts getBisulfiteCounts() {
        return bisulfiteCounts;
    }

    @Override
    public BaseModificationCounts getModifiedBaseCounts() {
        return baseModificationCounts;
    }

    @Override
    public int getBucketSize() {
        return bucketSize;
    }

    @Override
    public boolean hasBaseCounts() {
        return true;
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
        if (baseModificationCounts != null) {
            baseModificationCounts.incrementCounts(alignment);
        }

        int alignmentStart = alignment.getAlignmentStart();
        int alignmentEnd = alignment.getAlignmentEnd();
        Strand strand = alignment.getReadStrand();
        final boolean isNegativeStrand = strand == Strand.NEGATIVE;

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        if (blocks != null) {
            int lastBlockEnd = -1;

            for (AlignmentBlock b : blocks) {
                if (b.getEnd() < start) continue;
                if (b.getStart() > end) break;

                //Strand strand = alignment.getFirstOfPairStrand();

                // Don't count softclips
                if (!b.isSoftClip() && strand != Strand.NONE) {
                    incBlockCounts(b, isNegativeStrand);
                }
            }

            // Count deletions
            List<Gap> gaps = alignment.getGaps();
            if (gaps != null) {
                for (Gap gap : gaps) {
                    if (gap.getType() == SAMAlignment.DELETION) {
                        for (int pos = gap.getStart(); pos < gap.getStart() + gap.getnBases(); pos++) {
                            incrementDeletion(pos, isNegativeStrand);
                        }
                    }
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

        if (pos < getStart() || pos >= getEnd()) return null;

        StringBuffer buf = new StringBuffer();
        int totalCount = getTotalCount(pos);
        buf.append("Total count: " + totalCount);
        for (byte b : getBases()) {
            int negCount = getNegCount(pos, b);
            int posCount = getPosCount(pos, b);
            int count = negCount + posCount;
            int percent = Math.round(((float) count) * 100 / totalCount);
            char cU = Character.toUpperCase((char) b);
            buf.append("<br>" + cU + "      : " + count);
            if (count != 0) {
                buf.append("  (" + percent + "%,     " + posCount + "+,   " + negCount + "- )");
            }
        }
        int delCount = getDelCount(pos);
        int insCount = getInsCount(pos);
        if (delCount > 0 || insCount > 0) {
            buf.append("<br>DEL: " + delCount);
            buf.append("<br>INS: " + insCount);
        }

        return buf.toString();
    }


    /**
     * Return true if the mismatched (with respect to ref) read bases at the given position exceed the threshold.
     *
     * @param pos  genomic position (0 based)
     * @param ref  reference base
     * @param chr  chromosomes -- used as a key to fetch filtered snp locations
     * @param snpThreshold  threshold as a fraction of total
     * @return
     */
    public boolean isConsensusMismatch(int pos, byte ref, String chr, float snpThreshold) {

        boolean qualityWeight = PreferencesManager.getPreferences().getAsBoolean(SAM_ALLELE_USE_QUALITY);

        Set<Integer> filteredSnps = knownSnps == null ? null : knownSnps.get(chr);

        if (filteredSnps == null || !filteredSnps.contains(pos + 1)) {
            float threshold = snpThreshold * (qualityWeight ? getTotalQuality(pos) : getTotalCount(pos));
            float mismatchQualitySum = 0;

            if (ref > 0) {
                if (ref < 96) ref += 32;  // a fast "toLowercase"
                for (byte c : getBases()) {
                    if (c != ref && c != 'n') {
                        mismatchQualitySum += (qualityWeight ? getQuality(pos, (byte) c) : getCount(pos, (byte) c));
                    }
                }
                return (mismatchQualitySum >= threshold) && (threshold > 0); // (threshold > 0) avoids mismatch call in columns with all 0 quality
            }
        }
        return false;
    }

    public boolean isConsensusDeletion(int start, int width, float snpThreshold) {

        // We require deletion counts > threshold for at least 1/2 the width

        int end = start + width;
        int count = 0;
        for (int i = start; i < end; i++) {
            int totalCoverad = getTotalCount(i) + getDelCount(i);
            if (getDelCount(i) >= snpThreshold * totalCoverad) count++;
        }
        return count >= 0.5 * width;
    }

    @Override
    public boolean isConsensusInsertion(int pos, float snpThreshold) {
        float threshold = snpThreshold * (getTotalCount(pos) + getDelCount(pos)); // For this purpose consider deletions as covered
        return (this.getInsCount(pos) >= threshold);
    }


    /**
     * Load the set of known snps from a tab delimited file, format
     * chr < tab> location
     * The location is "1 base"  (first nucleotide is position 1).
     *
     * @param snpFile
     */
    private static void loadKnownSnps(String snpFile) {

        // This method might get called many times concurrently, but we only want to load these once.
        if (knownSnps != null) {
            return;
        }

        knownSnps = new HashMap<>();
        try (AsciiLineReader reader = ParsingUtils.openAsciiReader(new ResourceLocator(snpFile))) {
            String nextLine = "";
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0];
                Set<Integer> snps = knownSnps.computeIfAbsent(chr, k -> new HashSet<>(10000));
                snps.add(Integer.valueOf(tokens[1]));
            }
        } catch (Exception e) {
            knownSnps = null;
            log.error("", e);
            MessageUtils.showMessage("Error loading snps file: " + snpFile + " (" + e.toString() + ")");
        }


    }

    protected abstract void incPositionCount(int pos, byte n, byte q, boolean negativeStrand);

    protected abstract void incrementInsertion(AlignmentBlock insBlock);

    protected abstract void incrementDeletion(int pos, boolean negativeStrand);

    protected abstract void incBlockCounts(AlignmentBlock b, boolean isNegativeStrand);


}
