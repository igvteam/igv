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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

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

    public static final char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
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

    @Override
    public String getContig() {
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
                        if (gapTypes[gapIdx] == SAMAlignment.DELETION) {
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

        if (pos < getStart() || pos >= getEnd()) return null;

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
        int insCount = getInsCount(pos);
        buf.append("<br>---------------");
        if (delCount > 0 || insCount > 0) {
            buf.append("<br>DEL: " + delCount);
            buf.append("<br>INS: " + insCount);
        }

        return buf.toString();

    }

    public boolean isMismatch(int pos, byte ref, String chr, float snpThreshold) {

        boolean qualityWeight = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_ALLELE_USE_QUALITY);

        Set<Integer> filteredSnps = knownSnps == null ? null : knownSnps.get(chr);

        if (filteredSnps == null || !filteredSnps.contains(pos + 1)) {

            float threshold = snpThreshold * (qualityWeight ? getTotalQuality(pos) : getTotalCount(pos));
            float mismatchQualitySum = 0;

            if (ref > 0) {
                if (ref < 96) ref += 32;  // a fast "toLowercase"
                for (char c : nucleotides) {
                    if (c != ref && c != 'n') {
                        mismatchQualitySum += (qualityWeight ? getQuality(pos, (byte) c) : getCount(pos, (byte) c));
                    }

                }
                return mismatchQualitySum >= threshold;
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

    protected abstract void incBlockCounts(AlignmentBlock b, boolean isNegativeStrand);


}
