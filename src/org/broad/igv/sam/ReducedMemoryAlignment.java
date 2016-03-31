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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author jrobinso
 */
public class ReducedMemoryAlignment implements Alignment {

    private static Logger log = Logger.getLogger(ReducedMemoryAlignment.class);

    //String readName;
    private String chromosome;
    private int start;
    private int end;
    boolean negativeStrand;
    AlignmentBlock[] blocks;
    AlignmentBlock[] insertions;

    public ReducedMemoryAlignment(Alignment al, int indelLimit) {


        this.negativeStrand = al.isNegativeStrand();
        //  this.readName = al.getReadName();
        this.chromosome = al.getChr();
        this.start = al.getStart();
        this.end = al.getEnd();

        AlignmentBlock[] blocks = al.getAlignmentBlocks();
        if (blocks != null) {


            List<AlignmentBlock> rmBlocks = new ArrayList<AlignmentBlock>(blocks.length);
            int start = blocks[0].getStart();
            int end = blocks[0].getEnd();
            boolean softClip = blocks[0].isSoftClipped();

            for (int i = 1; i < blocks.length; i++) {

                if (blocks[i].getStart() - end < indelLimit && blocks[i].isSoftClipped() == softClip) {
                    end = blocks[i].getEnd();
                } else {
                    rmBlocks.add(new ReducedMemoryAlignmentBlock(start, end - start, softClip));
                    start = blocks[i].getStart();
                    end = blocks[i].getEnd();
                    softClip = blocks[i].isSoftClipped();
                }
            }

            // Last one
            rmBlocks.add(new ReducedMemoryAlignmentBlock(start, end - start, softClip));

            this.blocks = rmBlocks.toArray(new AlignmentBlock[rmBlocks.size()]);
        }

        AlignmentBlock[] insertions = al.getInsertions();
        if (insertions != null) {

            List<AlignmentBlock> rmInsertions = new ArrayList<AlignmentBlock>();
            for (AlignmentBlock b : insertions) {
                if (b.getLength() >= indelLimit) {
                    rmInsertions.add(b);
                }
            }
            this.insertions = rmInsertions.toArray(new AlignmentBlock[rmInsertions.size()]);
        }
    }


    public String getReadName() {
        return null;
    }

    /**
     * .aligned files do not include sequence
     *
     * @return
     */
    public String getReadSequence() {
        return "";
    }

    public void setMateSequence(String sequnce) {
        // Ignore
    }

    public String getPairOrientation() {
        return "";
    }

    public boolean isSmallInsert() {
        return false;
    }

    public boolean isVendorFailedRead() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public Color getColor() {
        return null;
    }


    public String getChromosome() {
        return chromosome;
    }

    public String getChr() {
        return chromosome;
    }

    @Override
    public String getContig() {
        return chromosome;
    }

    public int getAlignmentStart() {
        return getStart();
    }

    public boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    public AlignmentBlock[] getAlignmentBlocks() {
        return blocks;
    }

    public AlignmentBlock[] getInsertions() {
        return insertions;
    }

    public String getCigarString() {
        return "*";
    }

    public int getInferredInsertSize() {
        return 0;
    }

    public int getMappingQuality() {
        return 255;
    }

    public ReadMate getMate() {
        return null;
    }

    public boolean isProperPair() {
        return true;
    }

    public boolean isMapped() {
        return true;
    }

    public boolean isPaired() {
        return false;
    }

    public boolean isNegativeStrand() {
        return negativeStrand;
    }

    public boolean isDuplicate() {
        return false;
    }

    public float getScore() {
        return 1.0f;
    }

    public LocusScore copy() {
        return this;
    }

    public String getClipboardString(double location) {
        return getValueString(location, null);
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        return "Read length = " + (getEnd() - getStart());
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    public int getAlignmentEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    public byte getBase(double position) {
        return 0;
    }

    public byte getPhred(double position) {
        return 0;
    }

    public String getSample() {
        return null;
    }

    public String getReadGroup() {
        return null;
    }

    public String getLibrary() {
        return null;
    }

    public Object getAttribute(String key) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public List<Gap> getGaps() {
        return null;
    }

    public boolean isFirstOfPair() {
        return false;
    }

    public boolean isSecondOfPair() {
        return false;
    }

    public Strand getFirstOfPairStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }


    public Strand getSecondOfPairStrand() {
        return Strand.NONE;
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public void finish() {
    }

    @Override
    public boolean isPrimary() {
        return true;
    }

    @Override
    public boolean isSupplementary() {
        return false;
    }

    public static class ReducedMemoryAlignmentBlock implements AlignmentBlock {

        ReducedMemoryAlignmentBlock(int start, int length, boolean softClipped) {
            this.start = start;
            this.length = length;
            this.softClipped = softClipped;
        }

        int start;
        int length;
        boolean softClipped;

        @Override
        public boolean contains(int position) {
            int offset = position - start;
            return offset >= 0 && offset < getLength();
        }


        @Override
        public int getLength() {
            return length;
        }

        @Override
        public byte getBase(int offset) {
            return 0;
        }

        @Override
        public byte[] getBases() {
            return null;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public byte getQuality(int offset) {
            return 0;
        }

        @Override
        public byte[] getQualities() {
            return null;
        }

        @Override
        public int getEnd() {
            return start + length;
        }

        @Override
        public boolean isSoftClipped() {
            return softClipped;
        }

        @Override
        public void reduce(Genome genome) {

        }

        @Override
        public boolean hasBases() {
            return false;
        }

        @Override
        public FlowSignalSubContext getFlowSignalSubContext(int offset) {
            return null;
        }

        @Override
        public boolean hasFlowSignals() {
            return false;
        }
    }


    // TODO -- why is this class a "Feature"?

    public static class ReducedMemoryAlignmentCounts implements AlignmentCounts {

        private int start;
        private int end;
        private int bucketSize;
        private int nBuckets;
        double[] total;
        private int maxCount;


        public ReducedMemoryAlignmentCounts(int start, int end, int bucketSize) {
            this.start = start;
            this.bucketSize = bucketSize;

            this.nBuckets = (end - start) / bucketSize;
            if (nBuckets == 0 || (end - start) % bucketSize != 0) nBuckets++;

            total = new double[nBuckets];
            Arrays.fill(total, 0);

            // Adjust end to avoid fractional bucket
            this.end = start + nBuckets * bucketSize;

        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }

        @Override
        public int getBucketSize() {
            return bucketSize;
        }


        @Override
        public int getNumberOfPoints() {
            return nBuckets;
        }


        @Override
        public int getTotalCount(int pos) {
            int offset = (pos - start) / bucketSize;
            if (offset < 0 || offset >= total.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                return (int) Math.round(total[offset]);

            }
        }

        @Override
        public void incCounts(Alignment alignment) {

            AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
            if (blocks != null) {
                for (AlignmentBlock b : blocks) {
                    if(!b.isSoftClipped()) {
                        incrementBuckets(b.getStart(), b.getEnd());
                    }
                }

            } else {
                incrementBuckets(alignment.getAlignmentStart(), alignment.getAlignmentEnd());
            }
        }

        private void incrementBuckets(int blockStart, int blockEnd) {

            int startBucket = Math.max(0, (blockStart - this.start) / bucketSize);
            int endBucket = Math.min(nBuckets-1, (blockEnd - this.start) / bucketSize);

            for (int b = startBucket; b <= endBucket; b++) {
                int bucketStart = this.start + b * bucketSize;
                int bucketEnd = bucketStart + bucketSize;
                if(blockEnd >= bucketStart && blockStart <= bucketEnd) {
                    double s = Math.max(blockStart, bucketStart);
                    double e = Math.min(blockEnd, bucketEnd);
                    double f = (e - s) / bucketSize;
                    total[b] += f;

                    if(total[b] > maxCount) maxCount = (int) Math.round(total[b]);
                }

            }

        }


        @Override
        public int getMaxCount(int origin, int end) {
            return maxCount;
        }

        @Override
        public String getValueStringAt(int pos) {
            int idx = (pos - start) / bucketSize;
            return idx > 0 && idx < total.length ? String.valueOf((int) Math.round(total[idx])) : "";
        }

        // Rest is needed for the interface, but NA for reduced memory alignments
        @Override
        public String getChr() {
            return null;
        }

        @Override
        public String getContig() {
            return null;
        }


        @Override
        public int getTotalQuality(int pos) {
            return 0;
        }

        @Override
        public int getCount(int pos, byte b) {
            return 0;
        }

        @Override
        public int getNegCount(int pos, byte b) {
            return 0;
        }

        @Override
        public int getPosCount(int pos, byte b) {
            return 0;
        }

        @Override
        public int getDelCount(int pos) {
            return 0;
        }

        @Override
        public int getInsCount(int pos) {
            return 0;
        }

        @Override
        public int getQuality(int pos, byte b) {
            return 0;
        }


        @Override
        public boolean isMismatch(int pos, byte ref, String chr, float snpThreshold) {
            return false;
        }

        @Override
        public BisulfiteCounts getBisulfiteCounts() {
            return null;
        }

        @Override
        public boolean hasBaseCounts() {
            return false;
        }

        @Override
        public void finish() {

        }

    }
}
