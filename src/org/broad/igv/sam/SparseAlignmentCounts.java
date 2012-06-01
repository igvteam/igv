package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.util.collections.IntArrayList;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 11/22/11
 */
public class SparseAlignmentCounts extends BaseAlignmentCounts {

    private static Logger log = Logger.getLogger(SparseAlignmentCounts.class);
    private int maxCount = 0;
    List<Integer> indeces;


    /**
     * Map of genomic position -> index of count arrays
     */
    HashMap<Integer, Integer> indexMap;

    IntArrayList posA;
    IntArrayList posT;
    IntArrayList posC;
    IntArrayList posG;
    IntArrayList posN;
    IntArrayList negA;
    IntArrayList negT;
    IntArrayList negC;
    IntArrayList negG;
    IntArrayList negN;
    IntArrayList qA;
    IntArrayList qT;
    IntArrayList qC;
    IntArrayList qG;
    IntArrayList qN;
    IntArrayList posTotal;
    IntArrayList negTotal;
    IntArrayList del;
    IntArrayList ins;
    private IntArrayList totalQ;


    public SparseAlignmentCounts(int start, int end, AlignmentTrack.BisulfiteContext bisulfiteContext) {

        super(start, end, bisulfiteContext);

        indexMap = new HashMap<Integer, Integer>(1000);
        posA = new IntArrayList(1000);
        posT = new IntArrayList(1000);
        posC = new IntArrayList(1000);
        posG = new IntArrayList(1000);
        posN = new IntArrayList(1000);
        posTotal = new IntArrayList(1000);
        negA = new IntArrayList(1000);
        negT = new IntArrayList(1000);
        negC = new IntArrayList(1000);
        negG = new IntArrayList(1000);
        negN = new IntArrayList(1000);
        negTotal = new IntArrayList(1000);
        qA = new IntArrayList(1000);
        qT = new IntArrayList(1000);
        qC = new IntArrayList(1000);
        qG = new IntArrayList(1000);
        qN = new IntArrayList(1000);
        del = new IntArrayList(1000);
        ins = new IntArrayList(1000);
        totalQ = new IntArrayList(1000);
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getNumberOfPoints() {
        return indeces == null ? 0 : indeces.size();
    }

    public int getPosition(int idx) {
        return indeces.get(idx);
    }


    public int getMaxCount() {
        return maxCount;
    }

    public int getTotalCount(int pos) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            return getCountFromList(posTotal, idx) + getCountFromList(negTotal, idx);

        }
    }

    public int getNegTotal(int pos) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            return getCountFromList(negTotal, idx);

        }
    }

    public int getPosTotal(int pos) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            return getCountFromList(posTotal, idx);

        }
    }

    public int getTotalQuality(int pos) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            return getCountFromList(totalQ, idx);

        }
    }

    public int getCount(int pos, byte b) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            switch (b) {
                case 'a':
                case 'A':
                    return getCountFromList(posA, idx) + getCountFromList(negA, idx);
                case 't':
                case 'T':
                    return getCountFromList(posT, idx) + getCountFromList(negT, idx);
                case 'c':
                case 'C':
                    return getCountFromList(posC, idx) + getCountFromList(negC, idx);
                case 'g':
                case 'G':
                    return getCountFromList(posG, idx) + getCountFromList(negG, idx);
                case 'n':
                case 'N':
                    return getCountFromList(posN, idx) + getCountFromList(negN, idx);
            }
            log.debug("Unknown nucleotide: " + b);
            return 0;
        }
    }

    private int getCountFromList(IntArrayList list, int idx) {
        return idx < list.size() ? list.get(idx) : 0;
    }

    public int getNegCount(int pos, byte b) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            switch (b) {
                case 'a':
                case 'A':
                    return getCountFromList(negA, idx);
                case 't':
                case 'T':
                    return getCountFromList(negT, idx);
                case 'c':
                case 'C':
                    return getCountFromList(negC, idx);
                case 'g':
                case 'G':
                    return getCountFromList(negG, idx);
                case 'n':
                case 'N':
                    return getCountFromList(negN, idx);
            }
            log.error("Unknown nucleotide: " + b);
            return 0;
        }
    }

    public int getPosCount(int pos, byte b) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            switch (b) {
                case 'a':
                case 'A':
                    return getCountFromList(posA, idx);
                case 't':
                case 'T':
                    return getCountFromList(posT, idx);
                case 'c':
                case 'C':
                    return getCountFromList(posC, idx);
                case 'g':
                case 'G':
                    return getCountFromList(posG, idx);
                case 'n':
                case 'N':
                    return getCountFromList(posN, idx);
            }
            log.error("Unknown nucleotide: " + b);
            return 0;
        }
    }

    public int getDelCount(int pos) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            return getCountFromList(del, idx);
        }
    }


    public int getInsCount(int pos) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            return getCountFromList(ins, idx);
        }
    }

    public int getQuality(int pos, byte b) {
        if (!indexMap.containsKey(pos)) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            int idx = getIndex(pos);
            switch (b) {
                case 'a':
                case 'A':
                    return getCountFromList(qA, idx);
                case 't':
                case 'T':
                    return getCountFromList(qT, idx);
                case 'c':
                case 'C':
                    return getCountFromList(qC, idx);
                case 'g':
                case 'G':
                    return getCountFromList(qG, idx);
                case 'n':
                case 'N':
                    return getCountFromList(qN, idx);
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
        int idx = getIndex(pos);
        increment(del, idx, 1);
        if (countDeletedBasesCovered) {
            if (negativeStrand) {
                increment(negTotal, idx, 1);
            } else {
                increment(posTotal, idx, 1);
            }
        }
    }

    protected void incrementInsertion(AlignmentBlock insBlock) {
        int pos = insBlock.getStart();
        int idx1 = getIndex(pos);
        // Insertions are between bases.  increment count on either side
        increment(ins, idx1, 1);
        if (pos > 0) {
            int idx2 = getIndex(pos - 1);
            increment(ins, idx2, 1);
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

        int idx = getIndex(pos);
        switch (b) {
            case 'a':
            case 'A':
                if (isNegativeStrand) {
                    increment(negA, idx, 1);
                } else {
                    increment(posA, idx, 1);
                }
                increment(qA, idx, q);
                break;
            case 't':
            case 'T':
                if (isNegativeStrand) {
                    increment(negT, idx, 1);
                } else {
                    increment(posT, idx, 1);
                }
                increment(qT, idx, q);
                break;
            case 'c':
            case 'C':
                if (isNegativeStrand) {
                    increment(negC, idx, 1);
                } else {
                    increment(posC, idx, 1);
                }
                increment(qC, idx, q);
                break;
            case 'g':
            case 'G':
                if (isNegativeStrand) {
                    increment(negG, idx, 1);
                } else {
                    increment(posG, idx, 1);
                }
                increment(qG, idx, q);
                break;
            // Everything else is counted as "N".  This might be an actual "N",  or an ambiguity code
            default:
                if (isNegativeStrand) {
                    increment(negN, idx, 1);
                } else {
                    increment(posN, idx, 1);
                }
                increment(qN, idx, q);

        }

        if (isNegativeStrand) {
            increment(negTotal, idx, 1);
        } else {
            increment(posTotal, idx, 1);
        }
        increment(totalQ, idx, q);

        int pt = idx < posTotal.size() ? posTotal.get(idx) : 0;
        int nt = idx < negTotal.size() ? negTotal.get(idx) : 0;
        maxCount = Math.max(pt + nt, maxCount);

    }

    private int getIndex(int pos) {
        Integer index = indexMap.get(pos);
        if (index == null) {
            index = new Integer(indexMap.size());
            indexMap.put(pos, index);
        }
        return index.intValue();
    }


    private void increment(IntArrayList list, int idx, int delta) {
        if (idx < list.size()) {
            list.set(idx, list.get(idx) + delta);
        } else {
            list.set(idx, delta);
        }
    }

    public void finish() {
        indeces = new ArrayList<Integer>(indexMap.keySet());
        Collections.sort(indeces);

    }

}
