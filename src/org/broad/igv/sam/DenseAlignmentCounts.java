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

/**
 * @author jrobinso
 * @date Feb 23, 2011
 */
public class DenseAlignmentCounts extends BaseAlignmentCounts {

    private static Logger log = Logger.getLogger(DenseAlignmentCounts.class);

    // counts
    int[] posA;
    int[] posT;
    int[] posC;
    int[] posG;
    int[] posN;
    int[] negA;
    int[] negT;
    int[] negC;
    int[] negG;
    int[] negN;
    int[] qA;
    int[] qT;
    int[] qC;
    int[] qG;
    int[] qN;
    int[] posTotal;
    int[] negTotal;
    int[] del;
    int[] ins;
    private int[] totalQ;
    private int maxCount = 0;

    /**
     * We store the maximum number of counts over intervals
     * For autoscaling, doesn't have to be super precise
     */
    protected static int MAX_COUNT_INTERVAL = 100;
    protected int[] maxCounts;

    public DenseAlignmentCounts(int start, int end, AlignmentTrack.BisulfiteContext bisulfiteContext) {
        super(start, end, bisulfiteContext);

        int nPts = end - start;
        posA = new int[nPts];
        posT = new int[nPts];
        posC = new int[nPts];
        posG = new int[nPts];
        posN = new int[nPts];
        posTotal = new int[nPts];
        negA = new int[nPts];
        negT = new int[nPts];
        negC = new int[nPts];
        negG = new int[nPts];
        negN = new int[nPts];
        negTotal = new int[nPts];
        qA = new int[nPts];
        qT = new int[nPts];
        qC = new int[nPts];
        qG = new int[nPts];
        qN = new int[nPts];
        del = new int[nPts];
        ins = new int[nPts];
        totalQ = new int[nPts];

        maxCounts = new int[(nPts / MAX_COUNT_INTERVAL) + 1];
        log.debug("nPts: " + nPts + " maxCounts.length: " + maxCounts.length);
    }

    public int getNumberOfPoints() {
        return end - start;
    }

    @Override
    public int getMaxCount(int strt, int end) {

        if(maxCounts == null || maxCounts.length == 0) return 1;

        strt = Math.max(0, strt);
        end = Math.min(getEnd(), end);
        int startMCI = Math.max(0, (strt-this.start) / MAX_COUNT_INTERVAL);
        int endMCI = Math.max(0, (end-this.start) / MAX_COUNT_INTERVAL);
        endMCI = Math.min(endMCI, maxCounts.length - 1);
        int max = maxCounts[startMCI];
        for(int mci= startMCI+1; mci <= endMCI; mci++){
            max = Math.max(max, maxCounts[mci]);
        }
        return max;
    }

    public void finish() {
        // Noop
    }

    public int getTotalCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return posTotal[offset] + negTotal[offset];

        }
    }

    public int getTotalQuality(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return totalQ[offset];

        }
    }

    public int getCount(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            switch (b) {
                case 'a':
                case 'A':
                    return posA[offset] + negA[offset];
                case 't':
                case 'T':
                    return posT[offset] + negT[offset];
                case 'c':
                case 'C':
                    return posC[offset] + negC[offset];
                case 'g':
                case 'G':
                    return posG[offset] + negG[offset];
                case 'n':
                case 'N':
                    return posN[offset] + negN[offset];
            }
            log.debug("Unknown nucleotide: " + b);
            return 0;
        }
    }

    public int getNegCount(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            switch (b) {
                case 'a':
                case 'A':
                    return negA[offset];
                case 't':
                case 'T':
                    return negT[offset];
                case 'c':
                case 'C':
                    return negC[offset];
                case 'g':
                case 'G':
                    return negG[offset];
                case 'n':
                case 'N':
                    return negN[offset];
            }
            log.error("Unknown nucleotide: " + b);
            return 0;
        }
    }

    public int getPosCount(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            switch (b) {
                case 'a':
                case 'A':
                    return posA[offset];
                case 't':
                case 'T':
                    return posT[offset];
                case 'c':
                case 'C':
                    return posC[offset];
                case 'g':
                case 'G':
                    return posG[offset];
                case 'n':
                case 'N':
                    return posN[offset];
            }
            log.error("Unknown nucleotide: " + b);
            return 0;
        }
    }

    public int getDelCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        }
        return del[offset];
    }


    public int getInsCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        }
        return ins[offset];
    }

    public int getQuality(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            log.error("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            return 0;
        } else {
            switch (b) {
                case 'a':
                case 'A':
                    return qA[offset];
                case 't':
                case 'T':
                    return qT[offset];
                case 'c':
                case 'C':
                    return qC[offset];
                case 'g':
                case 'G':
                    return qG[offset];
                case 'n':
                case 'N':
                    return qN[offset];
            }
            log.error("Unknown nucleotide: " + posN);
            return 0;
        }
    }


    protected void incrementDeletion(int pos, boolean negativeStrand) {
        int offset = pos - start;
        if (offset >= 0 && offset < del.length) {
            del[offset] = del[offset] + 1;

            if (countDeletedBasesCovered) {
                if (negativeStrand) {
                    negTotal[offset] = negTotal[offset] + 1;
                } else {
                    posTotal[offset] = posTotal[offset] + 1;
                }
            }
        }
    }

    protected void incrementInsertion(AlignmentBlock insBlock) {
        int pos = insBlock.getStart();
        int offset = pos - start;
        // Insertions are between bases.  increment count on either side
        if (offset >= 0 && offset < ins.length) {
            ins[offset] = ins[offset] + 1;
            offset--;
            if (offset >= 0) {
                ins[offset] = ins[offset] + 1;
            }
        }
    }


    protected void incBlockCounts(AlignmentBlock block, boolean isNegativeStrand) {
        int start = block.getStart();
        byte[] bases = block.getBases();
        if (bases != null) {
            for (int i = 0; i < bases.length; i++) {
                int pos = start + i;
                // NOTE:  the direct access block.qualities is intentional,  profiling reveals this to be a critical bottleneck
                byte q = block.qualities[i];
                // TODO -- handle "=" in cigar string with no read bases
                byte n = bases[i];
                incPositionCount(pos, n, q, isNegativeStrand);
            }
        }
    }

    protected void incPositionCount(int pos, byte b, byte q, boolean isNegativeStrand) {

        int offset = pos - start;
        if (offset >= 0 && offset < posA.length) {
            switch (b) {
                case 'a':
                case 'A':
                    if (isNegativeStrand) {
                        negA[offset] = negA[offset] + 1;
                    } else {
                        posA[offset] = posA[offset] + 1;
                    }
                    qA[offset] = qA[offset] + q;
                    break;
                case 't':
                case 'T':
                    if (isNegativeStrand) {
                        negT[offset] = negT[offset] + 1;
                    } else {
                        posT[offset] = posT[offset] + 1;
                    }
                    qT[offset] = qT[offset] + q;
                    break;
                case 'c':
                case 'C':
                    if (isNegativeStrand) {
                        negC[offset] = negC[offset] + 1;
                    } else {
                        posC[offset] = posC[offset] + 1;
                    }
                    qC[offset] = qC[offset] + q;
                    break;
                case 'g':
                case 'G':
                    if (isNegativeStrand) {
                        negG[offset] = negG[offset] + 1;
                    } else {
                        posG[offset] = posG[offset] + 1;
                    }
                    qG[offset] = qG[offset] + q;
                    break;
                // Everything else is counted as "N".  This might be an actual "N",  or an ambiguity code
                default:
                    if (isNegativeStrand) {
                        negN[offset] = negN[offset] + 1;
                    } else {
                        posN[offset] = posN[offset] + 1;
                    }
                    qN[offset] = qN[offset] + q;

            }
            if (isNegativeStrand) {
                negTotal[offset] = negTotal[offset] + 1;
            } else {
                posTotal[offset] = posTotal[offset] + 1;
            }
            totalQ[offset] = totalQ[offset] + q;

            int tmp = posTotal[offset] + negTotal[offset];
            int maxCountInt = offset / MAX_COUNT_INTERVAL;
            if(tmp > maxCounts[maxCountInt]){
                maxCounts[maxCountInt] = tmp;
            }
        }
    }

}
