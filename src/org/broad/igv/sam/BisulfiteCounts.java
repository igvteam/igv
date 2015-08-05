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

import org.broad.igv.feature.genome.Genome;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 2/6/12
 */
public class BisulfiteCounts {

    public static final AlignmentTrack.BisulfiteContext[] NOMESEQ_CONTEXTS =
            {AlignmentTrack.BisulfiteContext.HCG, AlignmentTrack.BisulfiteContext.GCH};

    Genome genome;

    AlignmentTrack.BisulfiteContext bisulfiteContext;

    Map<Integer, Count> counts = new HashMap<Integer, Count>();

    public BisulfiteCounts(AlignmentTrack.BisulfiteContext bisulfiteContext, Genome genome) {
        this.bisulfiteContext = bisulfiteContext;
        this.genome = genome;
    }

    public void incrementCounts(Alignment baseAlignment) {

        // Only works with block formats
        if(baseAlignment.getAlignmentBlocks() == null) return;

        String chrname = genome.getChromosomeAlias(baseAlignment.getChr());

        boolean flipRead;
        // We will only need reverse complement if the strand and paired end status don't match (2nd ends are G->A)
        if (baseAlignment.isPaired()) {
            flipRead = (baseAlignment.isPaired() && (baseAlignment.isNegativeStrand() ^ baseAlignment.isSecondOfPair()));
        } else {
            flipRead = baseAlignment.isNegativeStrand();
        }


        for (AlignmentBlock block : baseAlignment.getAlignmentBlocks()) {

            if (block.isSoftClipped()) continue;

            int start = block.getStart();
            int end = block.getEnd();
            byte[] inReference = genome.getSequence(chrname, start, end);
            if (inReference == null) continue;
            byte[] reference = (flipRead ?
                    AlignmentUtils.reverseComplementCopy(inReference) :
                    inReference);

            byte[] inRead = block.getBases();
            byte[] read = (flipRead) ? AlignmentUtils.reverseComplementCopy(inRead) : inRead;

            int alignmentLen = inRead.length;
            final int idxEnd = alignmentLen - 1;
            for (int idxFw = 0; idxFw < alignmentLen; idxFw++) {

                //// Everything comes in relative to the positive strand.
                int idx = flipRead ? (idxEnd - idxFw) : idxFw;


                // Since we allow soft-clipping, the reference sequence can actually be shorter than the read.  Not sure
                // what to do in this case,  just skip?
                if (idx < 0 || idx >= reference.length) continue;

                // The read base can be an equals sign, so change that to the actual ref base
                byte refbase =  reference[idx];


                // Strand has already been accounted for
                if (refbase == 'C') {
                    byte readbase = read[idx];
                    if (readbase == '=') readbase = refbase;

                    if (AlignmentUtils.compareBases((byte) 'C', readbase) || AlignmentUtils.compareBases((byte) 'T', readbase)) {

                        AlignmentTrack.BisulfiteContext matchingContext = contextIsMatching(reference, read, idx, bisulfiteContext);
                        boolean matchesContext = (matchingContext != null);
                        if (matchesContext) {
                            //int extension = this.getBisulfiteSymmetricCytosineExtension(bisulfiteContext);
                            int extension = 0; // Don't extend
                            int step = flipRead ? -1 : 1;
                            for (int i = 0; i <= extension; i++) {
                                int position = start + idxFw + step*i;

                                Count count = counts.get(position);
                                if (count == null) {
                                    count = new Count();
                                    counts.put(position, count);
                                }
                                if (AlignmentUtils.compareBases((byte) 'T', readbase)) {
                                    count.increment(false);
                                } else if (AlignmentUtils.compareBases((byte) 'C', readbase)) {
                                    count.increment(true);
                                }
                            }
                        }
                    }
                }
                else {
                    // Non-informative
                }
            }
        }
    }

    public Count getCount(int position) {
        Count c = counts.get(position);
        if (c == null) {
            c = new Count();
            counts.put(position, c);
        }
        return c;
    }

    /**
     * @param reference
     * @param read
     * @param idx
     * @param bisulfiteContext
     * @return Returns the context that matched (in the case of the base class, this is always the same
     *         as the context passed in, derived classes might return a different context).
     *         If we don't match, return null.
     */
    protected AlignmentTrack.BisulfiteContext contextIsMatching(byte[] reference, byte[] read, int idx,
                                                                AlignmentTrack.BisulfiteContext bisulfiteContext) {


        // TODO -- NOMESEQ
//        for (AlignmentTrack.BisulfiteContext context : NOMESEQ_CONTEXTS) {
//            if (super.contextIsMatching(reference, read, idx, context) != null) return context;
//        }
//        return null;


        // Get the context and see if it matches our desired context.
        byte[] preContext = AlignmentTrack.getBisulfiteContextPreContext(bisulfiteContext);
        byte[] postContext = AlignmentTrack.getBisulfiteContextPostContext(bisulfiteContext);

        boolean matchesContext = true;

        // First do the "post" context
        int minLen = Math.min(reference.length, read.length);
        if ((idx + postContext.length) >= minLen) {
            matchesContext = false;
        } else {
            // Cut short whenever we don't match
            for (int posti = 0; matchesContext && (posti < postContext.length); posti++) {
                byte contextb = postContext[posti];
                int offsetidx = idx + 1 + posti;
                matchesContext &= positionMatchesContext(contextb, reference[offsetidx], read[offsetidx]);
            }
        }

        // Now do the pre context
        if ((idx - preContext.length) < 0) {
            matchesContext = false;
        } else {
            // Cut short whenever we don't match
            for (int prei = 0; matchesContext && (prei < preContext.length); prei++) {
                byte contextb = preContext[prei];
                int offsetidx = idx - (preContext.length - prei);
                matchesContext &= positionMatchesContext(contextb, reference[offsetidx], read[offsetidx]);
            }
        }

        return (matchesContext) ? bisulfiteContext : null;
    }

    /**
     * @param contextb      The residue in the context string (IUPAC)
     * @param referenceBase The reference sequence (already checked that offsetidx is within bounds)
     * @param readBase      The read sequence (already checked that offsetidx is within bounds)
     * @return
     */
    protected boolean positionMatchesContext(byte contextb, final byte referenceBase, final byte readBase) {

        boolean matchesContext = AlignmentUtils.compareBases(contextb, referenceBase);
        if (!matchesContext) {
            return false; // Don't need to check any further
        }

        // For the read, we have to handle C separately
        boolean matchesReadContext = AlignmentUtils.compareBases(contextb, readBase);
        if (AlignmentUtils.compareBases((byte) 'T', readBase)) {
            matchesReadContext |= AlignmentUtils.compareBases(contextb, (byte) 'C');
        }

        return matchesReadContext;
    }


    /**
     * caller must be careful to shift relative to the strand the cytosine is on).
     */
    protected int getBisulfiteSymmetricCytosineExtension(AlignmentTrack.BisulfiteContext item) {

        int out;
        switch (item) {
            case CG:
            case HCG:
            case WCG:
                out = 1;
                break;
            case CHG:
                out = 2;
                break;
            case GCH:   // Added by JTR,  confirm?
                out = -1;
                break;
            default:
                out = 0;
                break;
        }

        return out;
    }


    public static class Count {
        int methylatedCount;
        int unmethylatedCount;

        public void increment(boolean methylated) {
            if (methylated) methylatedCount++;
            else unmethylatedCount++;
        }
    }

}
