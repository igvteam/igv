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

package org.broad.igv.maf;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 2/18/13
 *         Time: 9:49 AM
 */
public class MultipleAlignmentBlock {

    private String chr;
    private int start;
    private int end;
    private int[] gapAdjustedIndex;
    private double score;
    private List<Sequence> sequences;
    private List<Gap> gaps;
    private String key;

    public MultipleAlignmentBlock() {
        sequences = new ArrayList<Sequence>();
        gaps = new ArrayList<Gap>();
    }


    /**
     * Return a key to uniquely identify this alignment.
     *
     * @return
     */
    public String getKey() {
        if (key == null) {
            StringBuilder b = new StringBuilder();
            b.append(chr);
            b.append(":");
            b.append(start);
            b.append("-");
            b.append(end);
            for (Sequence seq : sequences) {
                b.append("_");
                b.append(seq.getSpecies());
            }
            key = b.toString();
        }
        return key;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public double getScore() {
        return score;
    }

    public List<Sequence> getSequences() {
        return sequences;
    }

    public int getGapAdjustedIndex(int position) {
        return gapAdjustedIndex[position - start];
    }


    public Sequence getRefSequence() {
        return sequences.get(0);
    }

    public List<Gap> getGaps() {
        return gaps;
    }

    public void addSequence(Sequence sequence) {
        if (sequences.isEmpty()) {
            // First sequence is the reference
            chr = sequence.chr;
            start = sequence.start;
            end = start + sequence.size;
            findGaps(sequence.text);
        }
        sequences.add(sequence);
    }

    private void findGaps(String referenceText) {

        byte[] bytes = referenceText.getBytes();

        // Maps genomic position -> text position.
        gapAdjustedIndex = new int[end - start];
        int refPosition = 0;

        byte dash = (byte) '-';

        Gap gap = null;

        for (int i = 0; i < bytes.length; i++) {
            if (bytes[i] == dash) {
                if (gap == null) {
                    // Start a new gap.
                    double gapPosition = start + refPosition;
                    gap = new Gap(gapPosition, i);
                } else {
                    // Increment existing gap
                    gap.incSize();
                }
            } else {

                if (gap != null) {
                    // End in-progress gap
                    gaps.add(gap);
                    gap = null;
                }
                gapAdjustedIndex[refPosition] = i;
                refPosition++;
            }
        }

        if (gap != null) {
            gaps.add(gap);
        }


    }

    public Sequence getSequence(String sp) {
        // TODO -- use a map
        for(Sequence seq : sequences) {
            if(seq.getSpecies().equals(sp)) return seq;
        }
        return null;
    }


    /**
     * Class to represent a gap in the reference sequence.
     */
    public static class Gap {
        double position;
        int startIdx;
        int size;

        public Gap(double position, int startIdx) {
            this.position = position;
            this.startIdx = startIdx;
            this.size = 1;
        }

        public void incSize() {
            size++;
        }
    }

    public static class Sequence {

        private String species;
        private String chr;
        private int start;
        private int size;
        private char strand;
        private int srcSize;
        private String text;

        public Sequence(String species, String chr, int start, int size, char strand, int srcSize, String text) {
            this.species = species;
            this.chr = chr;
            this.start = start;
            this.size = size;
            this.strand = strand;
            this.srcSize = srcSize;
            this.text = text;
        }

        public String getSpecies() {
            return species;
        }

        public void setSpecies(String species) {
            this.species = species;
        }

        public String getChr() {
            return chr;
        }

        public void setChr(String chr) {
            this.chr = chr;
        }

        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public int getSize() {
            return size;
        }

        public void setSize(int size) {
            this.size = size;
        }

        public char getStrand() {
            return strand;
        }

        public void setStrand(char strand) {
            this.strand = strand;
        }

        public int getSrcSize() {
            return srcSize;
        }

        public void setSrcSize(int srcSize) {
            this.srcSize = srcSize;
        }

        public String getText() {
            return text;
        }

        public void setText(String text) {
            this.text = text;
        }
    }


}
