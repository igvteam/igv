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

import org.broad.igv.sam.AlignmentTrack.BisulfiteContext;

import java.awt.*;

/**
 * @author benb
 *         Benjamin Berman, University of Southern California
 *         <p/>
 *         Note that this is only currently supporting a single bisulfite protocol with the following assumptions:
 *         - The first end of paired end reads have C->T variants, while the second ends have G->A
 *         - Bisulfite only affects one strand.  So for the first end, you don't see G->A and for the
 *         second end you don't see C->T.  This allows us to detect true C->T genomic changes by examining
 *         the strand opposite the bisulfite event (we use the DEAMINATION_COLOR for these variants).
 *         <p/>
 *         This is the Illumina protocol published by Joe Ecker lab (Lister et al. 2009, Lister et al. 2011)
 *         and our lab (Berman et al. 2011)
 */
public class BisulfiteBaseInfo {


    public enum DisplayStatus {
        NOTHING, COLOR, CHARACTER
    }

    // Constants
    public static Color CYTOSINE_MISMATCH_COLOR = new Color(139, 94, 60); // Brown 
    public static Color NONCYTOSINE_MISMATCH_COLOR = new Color(139, 94, 60); // Brown 
    public static Color DEAMINATION_COLOR = new Color(139, 94, 60); // Brown 
    // Color.BLACK; // Like to make this brown but it doesn't scale with quality // 
    public static Color METHYLATED_COLOR = Color.red;
    public static Color UNMETHYLATED_COLOR = Color.blue;


    // Private vars
    private DisplayStatus[] displayStatus = null;
    private byte[] displayChars = null;
    private Color[] displayColors = null;
    private BisulfiteContext myContext = null;
    private boolean flipRead;


    /**
     * We require loading this once for the entire alignment so that we don't have to calculate reverse complements
     * more than once.
     *
     * @param inReference
     * @param block
     * @param bisulfiteContext
     */
    public BisulfiteBaseInfo(byte[] inReference, Alignment baseAlignment, AlignmentBlock block, BisulfiteContext bisulfiteContext) {
        super();

        myContext = bisulfiteContext;
        byte[] inRead = block.getBases();
        int alignmentLen = inRead.length;

        // We will only need reverse complement if the strand and paired end status don't match (2nd ends are G->A)
        if (baseAlignment.isPaired()) {
            flipRead = (baseAlignment.isPaired() && (baseAlignment.isNegativeStrand() ^ baseAlignment.isSecondOfPair()));
        } else {
            flipRead = baseAlignment.isNegativeStrand();
        }
        byte[] read = (flipRead) ? AlignmentUtils.reverseComplementCopy(inRead) : inRead;
        byte[] reference = (flipRead && inReference != null) ? AlignmentUtils.reverseComplementCopy(inReference) : inReference;


        displayChars = new byte[alignmentLen];
        displayStatus = new DisplayStatus[alignmentLen];
        displayColors = new Color[alignmentLen];

        final int idxEnd = alignmentLen - 1;
        for (int idxFw = 0; idxFw < alignmentLen; idxFw++) {

            //// Everything comes in relative to the positive strand.
            int idx = flipRead ? (idxEnd - idxFw) : idxFw;


            // Since we allow soft-cliping, the reference sequence can actually be shorter than the read.  Not sure
            // what to do in this case,  just skip?
            if (idx < 0 || idx > reference.length) continue;

            // The read base can be an equals sign, so change that to the actual ref base
            byte refbase = reference[idx];
            byte readbase = read[idx];
            if (readbase == '=') readbase = refbase;

            Color out = null;
            boolean matchesContext = false;


            // The logic is organized according to the reference base.

            switch (refbase) {
                case 'T':
                case 'A':
                case 'G':
                	if (AlignmentUtils.compareBases((byte) 'C', readbase)) out = METHYLATED_COLOR;
                	else if (!AlignmentUtils.compareBases(readbase, refbase)) out = NONCYTOSINE_MISMATCH_COLOR;
                    break;
                case 'C':
                    if (!AlignmentUtils.compareBases((byte) 'C', readbase) && !AlignmentUtils.compareBases((byte) 'T', readbase)) {
                        out = CYTOSINE_MISMATCH_COLOR;
                    } else {
                    	// If we had information about whether this position was a SNP or not, we could
                    	// show cytosines in any context when they are a SNP.
                        BisulfiteContext matchingContext = contextIsMatching(reference, read, idx, bisulfiteContext);
                        matchesContext = (matchingContext != null);
                        if (matchesContext) {
                            out = getContextColor(readbase, matchingContext);
                        }
                    }
                    break;
//                case 'G':
//                    if (AlignmentUtils.compareBases((byte) 'A', readbase))
//                    {
//                    	out = DEAMINATION_COLOR;
//                    }
//                    else if (!AlignmentUtils.compareBases(readbase, refbase))
//                    {
//                    	out = NONCYTOSINE_MISMATCH_COLOR;
//                    }
            }

            // Remember, the output should be relative to the FW strand (use idxFw)
            this.displayColors[idxFw] = out;
            if (out == null) {
                this.displayStatus[idxFw] = DisplayStatus.NOTHING;
            } else {
                if (matchesContext) {
                    // Display the color
                    this.displayStatus[idxFw] = DisplayStatus.COLOR;
                } else {
                    // Display the character
                    this.displayStatus[idxFw] = DisplayStatus.CHARACTER;
                    this.displayChars[idxFw] = 'X';
                }
            }
//			System.err.printf("\tSeting displayStatus[%d] = %s\n", idx, displayStatus[idx]);
        }

        //this.numDisplayStatus();
    }


    protected Color getContextColor(byte readbase,
                                    BisulfiteContext bisulfiteContext) {
        Color out = null;
        if (AlignmentUtils.compareBases((byte) 'T', readbase)) {
            out = UNMETHYLATED_COLOR;
        } else if (AlignmentUtils.compareBases((byte) 'C', readbase)) {
            out = METHYLATED_COLOR;
        }
        // C and T should be the only options at this point
        //						else
        //						{
        //							out = Color.yellow;
        //						}

        return out;
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
    protected BisulfiteContext contextIsMatching(byte[] reference, byte[] read, int idx,
                                                 BisulfiteContext bisulfiteContext) {

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

                //				System.err.printf("POST posMatchesContext(posti=%d, contextb=%c, refb=%c, readb=%c, offsetidx=%d) = %s\n",
                //						posti, contextb, reference[offsetidx], read[offsetidx], offsetidx, matchesContext);

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
                //				System.err.printf("PRE posMatchesContext(prei=%d, contextb=%c, refb=%c, readb=%c, offsetidx=%d) = %s\n",
                //						prei, contextb, reference[offsetidx], read[offsetidx], offsetidx, matchesContext);
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

    public Color getDisplayColor(int idx) {
        return displayColors[idx];
    }

    public DisplayStatus getDisplayStatus(int idx) {
        return displayStatus[idx];
    }


    public int numDisplayStatus() {
        int len = this.displayStatus.length;
        boolean done = false;
        int i = 0;
        for (i = 0; (!done && (i < len)); i++) {
            done = (this.displayStatus[i] == null);
        }
//		System.err.printf("numDisplayStatus, len=%d\ti=%d\n", len, i);
        return i;
    }

    public double getXaxisShift(int idx) {

        double offset = 0.0;
        // This gets CpGs on opposite strands to line up. Only do it for meth values though, not snps

        if (getDisplayStatus(idx).equals(DisplayStatus.COLOR)) {
            double baseOffset = getBisulfiteSymmetricCytosineShift(myContext);
            offset = offset + ((flipRead ? -1 : 1) * baseOffset);
        }
        return offset;
    }

    /**
     * @param item
     * @return 0 if the cytosine context is not symmetric, else the number of
     *         base pairs to shift the position by (caller must be careful to shift
     *         relative to the strand the cytosine is on).
     */
    protected double getBisulfiteSymmetricCytosineShift(BisulfiteContext item) {
        double out = 0.0;

//        // The following may be too non-intuitive? BPB
//        switch (item) {
//            case CG:
//            case HCG:
//            case WCG:
//                out = 0.5;
//                break;
//            case CHG:
//                out = 1.0;
//                break;
//            case GCH:   // Added by JTR,  confirm?
//                out = -0.5;
//                break;
//            default:
//                out = 0.0;
//                break;
//        }

        return out;
    }


}
