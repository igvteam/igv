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

package org.broad.igv.synteny;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;


/**
 * // TODO -- consider implementing feature interface
 *
 * @author jrobinso
 */
public class BlastMapping extends BasicFeature {

    private Block queryBlock;
    private Block subjectBlock;

    /**
     * Constructs ...
     *
     * @param queryBlock
     * @param subjectBlock
     * @param percentIden
     */
    public BlastMapping(Block queryBlock, Block subjectBlock, float percentIden) {
        this.queryBlock = queryBlock;
        this.subjectBlock = subjectBlock;

        this.setStart(Math.min(queryBlock.getStart(), queryBlock.getEnd()));
        this.setEnd(Math.max(queryBlock.getStart(), queryBlock.getEnd()));
        this.setName("[" + subjectBlock.getContig() + ":" + subjectBlock.getStart() + "-" + subjectBlock.getEnd() + "]");
        this.setScore(percentIden);

    }


    public Block getQueryBlock() {
        return queryBlock;
    }


    public Block getSubjectBlock() {
        return subjectBlock;
    }

    @Override
    public String getChr() {
        return queryBlock.getContig();
    }

    public Strand getStrand() {
        return subjectBlock.getStrand();
    }

    public boolean containsQueryPosition(String contig, int position) {
        return queryBlock.containsPosition(contig, position);
    }


    /**
 * @author jrobinso
     */
    public static class Block {

        private String contig;
        private int start;
        private int end;
        Strand strand;

        /**
         * Constructs ...
         *
         * @param contig
         * @param start
         * @param end
         */
        public Block(String contig, int start, int end, Strand strand) {
            this.contig = contig;
            this.start = start;
            this.end = end;
            this.strand = strand;
        }


        public String getContig() {
            return contig;
        }


        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public Strand getStrand() {
            return strand;
        }

        public boolean containsPosition(String contig, int position) {
            return this.contig.equals(contig) && position >= start && position < end;
        }

    }
}
