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
    }

    public int getNumberOfPoints() {
        return end - start;
    }

    public void finish() {
        // Noop
    }

    public int getMaxCount() {
        return maxCount;
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

    public int getNegTotal(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return negTotal[offset];

        }
    }

    public int getPosTotal(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= posA.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return posTotal[offset];

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

    public int getAvgQuality(int pos, byte b) {
        int count = getCount(pos, b);
        return count == 0 ? 0 : getQuality(pos, b) / count;
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
            maxCount = tmp > maxCount ? tmp : maxCount;
        }
    }

    public AlignmentCounts merge(AlignmentCounts other, AlignmentTrack.BisulfiteContext bisulfiteContext) {
        if (other.getClass() != this.getClass()) {
            throw new IllegalArgumentException("Cannot merge different types of AlignmentCount instances");
        }
        return DenseAlignmentCounts.merge(this, (DenseAlignmentCounts) other, bisulfiteContext);
    }

    /**
     * Merge first and second into a third DenseAlignmentCounts
     * Input arguments unchanged
     *
     * @param first
     * @param second
     * @param bisulfiteContext
     * @return
     */
    private static DenseAlignmentCounts merge(DenseAlignmentCounts first, DenseAlignmentCounts second, AlignmentTrack.BisulfiteContext bisulfiteContext) {
        if (second.getStart() < first.getStart()) {
            DenseAlignmentCounts tmp = first;
            first = second;
            second = tmp;
        }

        //It is now guaranteed that first.getStart() < second.getStart()
        //However, the relationship between first.getEnd() and second.getEnd() is unknown
        int firstLength = first.getEnd() - first.getStart();

        int newStart = first.getStart();
        int newEnd = Math.max(first.getEnd(), second.getEnd());

        int lengthIncrease = (newEnd - newStart) - (firstLength);

        //No need for copying if second is a subset of the first
        if (lengthIncrease <= 0) {
            return first;
        }

        int[][] firstSrcArrs = getCountArrs(first);
        int[][] secondSrcArrs = getCountArrs(second);

        DenseAlignmentCounts result = new DenseAlignmentCounts(newStart, newEnd, bisulfiteContext);
        int[][] destArrs = getCountArrs(result);
        int number = firstSrcArrs.length;
        int secondOffset = first.getEnd() - second.getStart();
        int[] destArr;
        for (int arnum = 0; arnum < number; arnum++) {
            destArr = destArrs[arnum];
            System.arraycopy(firstSrcArrs[arnum], 0, destArr, 0, firstLength);
            System.arraycopy(secondSrcArrs[arnum], secondOffset, destArr, firstLength, lengthIncrease);
        }
        result.maxCount = Math.max(first.getMaxCount(), second.getMaxCount());
        return result;
    }

    /**
     * Create an array of the arrays we use to keep track of data
     *
     * @param counts
     * @return
     */
    private static int[][] getCountArrs(DenseAlignmentCounts counts) {

        int[][] result = {counts.posA, counts.posT, counts.posC, counts.posG, counts.posN,
                counts.negA, counts.negT, counts.negC, counts.negG, counts.negN,
                counts.qA, counts.qT, counts.qC, counts.qG, counts.qN, counts.posTotal, counts.negTotal,
                counts.del, counts.ins, counts.totalQ};
        return result;
    }


}
