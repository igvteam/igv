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

import org.broad.igv.logging.*;

import java.util.*;

/**
 * @author jrobinso
 * @date Feb 23, 2011
 */
public class DenseAlignmentCounts extends BaseAlignmentCounts {

    private static Logger log = LogManager.getLogger(DenseAlignmentCounts.class);

    private final Set<Byte> nucleotides;
    private final Map<Byte, int[]> posCounts;
    private final Map<Byte, int[]> negCounts;
    private final Map<Byte, int[]> qualities;

    // counts
    int nPts;
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

        nPts = end - start;
        nucleotides = new LinkedHashSet<>(List.of((byte) 'A', (byte) 'T', (byte) 'C', (byte) 'G', (byte) 'N'));

        posCounts = new HashMap<>();
        negCounts = new HashMap<>();
        qualities = new HashMap<>();
        for (byte nt : nucleotides) {
            posCounts.put(nt, new int[nPts]);
            negCounts.put(nt, new int[nPts]);
            qualities.put(nt, new int[nPts]);
        }
        posTotal = new int[nPts];
        negTotal = new int[nPts];

        del = new int[nPts];
        ins = new int[nPts];
        totalQ = new int[nPts];

        maxCounts = new int[(nPts / MAX_COUNT_INTERVAL) + 1];
        log.debug("nPts: " + nPts + " maxCounts.length: " + maxCounts.length);
    }

    public int getNumberOfPoints() {
        return end - start;
    }

    public Set<Byte> getBases() {
        return nucleotides;
    }

    @Override
    public int getMaxCount(int strt, int end) {

        if (maxCounts == null || maxCounts.length == 0) return 1;

        strt = Math.max(0, strt);
        end = Math.min(getEnd(), end);
        int startMCI = Math.max(0, (strt - this.start) / MAX_COUNT_INTERVAL);
        int endMCI = Math.max(0, (end - this.start) / MAX_COUNT_INTERVAL);
        endMCI = Math.min(endMCI, maxCounts.length - 1);


        int max = 1;
        for (int mci = startMCI; mci <= endMCI; mci++) {
            if (mci >= maxCounts.length) {
                log.error("startMCI index out of range: " + mci + " startMCI=" + startMCI + "  endMCI=" + endMCI);
                return max;
            }
            max = Math.max(max, maxCounts[mci]);
        }
        return max;
    }

    public void finish() {
        // Noop
    }

    public int getTotalCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return posTotal[offset] + negTotal[offset];
        }
    }

    public int getTotalPositiveCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return posTotal[offset];
        }
    }

    public int getTotalNegativeCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return negTotal[offset];
        }
    }

    public int getTotalQuality(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
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
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            b = toUpperCase(b);
            if (!posCounts.containsKey(b) || !negCounts.containsKey(b)) {
                return 0;
            }
            return posCounts.get(b)[offset] + negCounts.get(b)[offset];
        }
    }

    public int getNegCount(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            b = toUpperCase(b);  //   "toUppercase"
            return negCounts.get(b)[offset];
        }
    }

    public int getPosCount(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            b = toUpperCase(b);
            return posCounts.get(b)[offset];
        }
    }

    public int getDelCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        }
        return del[offset];
    }


    public int getInsCount(int pos) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        }
        return ins[offset];
    }

    public int getQuality(int pos, byte b) {
        int offset = pos - start;
        if (offset < 0 || offset >= nPts) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 32;
        } else {
            b = toUpperCase(b);
            if (qualities.containsKey(b)) {
                return qualities.get(b)[offset];
            } else {
                log.error("Unknown nucleotide: " + b);
                return 0;
            }
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
        // Insertions are between bases.  increment count at position just before insertion
        if (offset >= 0 && offset < ins.length) {
            ins[offset] = ins[offset] + 1;
        }
    }


    protected void incBlockCounts(AlignmentBlock block, boolean isNegativeStrand) {
        int start = block.getStart();
        ByteSubarray bases = block.getBases();
        if (bases != null) {
            for (int i = 0; i < bases.length; i++) {
                int pos = start + i;
                byte q = block.getQuality(i);
                byte n = bases.getByte(i);
                incPositionCount(pos, n, q, isNegativeStrand);
            }
        }
    }

    protected void incPositionCount(int pos, byte b, byte q, boolean isNegativeStrand) {

        int offset = pos - start;
        if (offset >= 0 && offset < nPts) {
            b = toUpperCase(b);
            if (!nucleotides.contains(b)) {
                nucleotides.add(b);
                posCounts.put(b, new int[nPts]);
                negCounts.put(b, new int[nPts]);
                qualities.put(b, new int[nPts]);
            }
            if (isNegativeStrand) {
                negCounts.get(b)[offset] = negCounts.get(b)[offset] + 1;
                negTotal[offset] = negTotal[offset] + 1;
            } else {
                posCounts.get(b)[offset] = posCounts.get(b)[offset] + 1;
                posTotal[offset] = posTotal[offset] + 1;
            }
            qualities.get(b)[offset] = qualities.get(b)[offset] + q;

            totalQ[offset] = totalQ[offset] + q;

            int tmp = posTotal[offset] + negTotal[offset];
            int maxCountInt = offset / MAX_COUNT_INTERVAL;
            if (tmp > maxCounts[maxCountInt]) {
                maxCounts[maxCountInt] = tmp;
            }
        }
    }

    private static byte toUpperCase(byte b) {
        if (b >= 'a' && b <= 'z') {
            return (byte) (b - ('a' - 'A'));
        }
        return b;
    }
}


