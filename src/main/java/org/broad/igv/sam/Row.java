package org.broad.igv.sam;

import java.util.ArrayList;
import java.util.List;

/**
 * A row of alignments, packed to minimize empty space
 */
public class Row  {

    int nextIdx;
    private double score = 0;
    private int lastEnd = Integer.MIN_VALUE;  // Track rightmost position for quick checks

    public List<Alignment> getAlignments() {
        return alignments;
    }

    List<Alignment> alignments;
    public double y;
    public double h;

    public Row() {
        nextIdx = 0;
        this.alignments = new ArrayList(100);
    }

    public void addAlignment(Alignment alignment) {
//        if (alignment instanceof LinkedAlignment) {
//            for (Alignment a : ((LinkedAlignment) alignment).alignments) {
//                alignments.add(a);
//            }
//        } else {
        alignments.add(alignment);
        // Update lastEnd to track the rightmost position in this row
        lastEnd = Math.max(lastEnd, alignment.getEnd());
//        }
    }

    /**
     * Get the rightmost end position of alignments in this row
     */
    public int getLastEnd() {
        return lastEnd;
    }

    public Alignment nextAlignment() {
        if (nextIdx < alignments.size()) {
            Alignment tmp = alignments.get(nextIdx);
            nextIdx++;
            return tmp;
        } else {
            return null;
        }
    }

    public int getNextStartPos() {
        if (nextIdx < alignments.size()) {
            return alignments.get(nextIdx).getStart();
        } else {
            return Integer.MAX_VALUE;
        }
    }

    public boolean hasNext() {
        return nextIdx < alignments.size();
    }

    public void resetIdx() {
        nextIdx = 0;
    }

    public double getScore() {
        return score;
    }

    //        @Override
//        public boolean equals(Object object){
//            if(!(object instanceof Row)){
//                return false;
//            }
//            Row other = (Row) object;
//            boolean equals = this.getStart() == other.getStart();
//            equals &= this.getLastEnd() == other.getLastEnd();
//            equals &= this.getScore() == other.getScore();
//
//            return equals;
//
//        }
//
//        @Override
//        public int hashCode(){
//            int score = (int) getScore();
//            score = score != 0 ? score : 1;
//            return (getStart() * getLastEnd() * score);
//        }

}
