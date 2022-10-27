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

import java.util.ArrayList;
import java.util.List;

/**
 * A row of alignments, packed to minimize empty space
 */
public class Row  {

    int nextIdx;
    private double score = 0;

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
//        }
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
