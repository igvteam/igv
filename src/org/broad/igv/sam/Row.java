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

import org.broad.igv.feature.Strand;

import java.util.ArrayList;
import java.util.List;

/**
 * A row of alignments, packed to minimize empty space
 *
 * @author jacob
 * @date 2014-Jan-10
 */
public class Row implements Comparable<Row> {
    int nextIdx;
    private double score = 0;
    List<Alignment> alignments;

    public Row() {
        nextIdx = 0;
        this.alignments = new ArrayList(100);
    }

    public void addAlignment(Alignment alignment) {
        alignments.add(alignment);
    }

    public void updateScore(AlignmentTrack.SortOption option, double center, AlignmentInterval interval, String tag) {
        setScore(calculateScore(option, center, interval, tag));
    }


    public double calculateScore(AlignmentTrack.SortOption option, double center, AlignmentInterval interval, String tag) {

        int adjustedCenter = (int) center;
        Alignment centerAlignment = AlignmentInterval.getFeatureContaining(alignments, adjustedCenter);
        if (centerAlignment == null) {
            return Integer.MAX_VALUE;
        } else {
            switch (option) {
                case START:
                    return centerAlignment.getStart();
                case STRAND:
                    return centerAlignment.isNegativeStrand() ? -1 : 1;
                case FIRST_OF_PAIR_STRAND:
                    Strand strand = centerAlignment.getFirstOfPairStrand();
                    int score = 2;
                    if (strand != Strand.NONE) {
                        score = strand == Strand.NEGATIVE ? 1 : -1;
                    }
                    return score;
                case NUCLEOTIDE:
                    byte base = centerAlignment.getBase(adjustedCenter);
                    byte ref = interval.getReference(adjustedCenter);

                    // Check insertions
                    int insertionScore = 0;
                    AlignmentBlock[] insertions = centerAlignment.getInsertions();
                    for (AlignmentBlock ins : insertions) {
                        int s = ins.getStart();
                        if (s == adjustedCenter || (s - 1) == adjustedCenter) {
                            insertionScore += ins.getBases().length;
                        }
                    }

                    float baseScore;
                    if (base == 'N' || base == 'n') {
                        baseScore = 2;  // Base is "n"
                    } else if (base == ref) {
                        baseScore = 3;  // Base is reference
                    } else {
                        //If base is 0, base not covered (splice junction) or is deletion.
                        if (base == 0) {
                            int delCount = interval.getDelCount(adjustedCenter);
                            if (delCount > 0) {
                                baseScore = -delCount;
                            } else {
                                //Base not covered, NOT a deletion
                                baseScore = 1;
                            }
                        } else {
                            int count = interval.getCount(adjustedCenter, base);
                            byte phred = centerAlignment.getPhred(adjustedCenter);
                            baseScore = -(count + (phred / 1000.0f));   // The second bit will always be < 1
                        }


                    }

                    return baseScore - insertionScore;

                case QUALITY:
                    return -centerAlignment.getMappingQuality();
                case SAMPLE:
                    String sample = centerAlignment.getSample();
                    score = sample == null ? 0 : sample.hashCode();
                    return score;
                case READ_GROUP:
                    String readGroup = centerAlignment.getReadGroup();
                    score = readGroup == null ? 0 : readGroup.hashCode();
                    return score;
                case INSERT_SIZE:
                    return -Math.abs(centerAlignment.getInferredInsertSize());
                case MATE_CHR:
                    ReadMate mate = centerAlignment.getMate();
                    if (mate == null) {
                        return Integer.MAX_VALUE;
                    } else {
                        if (mate.getChr().equals(centerAlignment.getChr())) {
                            return Integer.MAX_VALUE - 1;
                        } else {
                            return mate.getChr().hashCode();
                        }
                    }
                case SUPPLEMENTARY:
                    return centerAlignment.isSupplementary() ? 0 : 1;
                case TAG:
                    Object tagValue = centerAlignment.getAttribute(tag);
                    score = tagValue == null ? 0 : tagValue.hashCode();
                    return score;
                default:
                    return Integer.MAX_VALUE;
            }
        }
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

    public void setScore(double score) {
        this.score = score;
    }

    @Override
    public int compareTo(Row o) {
        return (int) Math.signum(getScore() - o.getScore());
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
