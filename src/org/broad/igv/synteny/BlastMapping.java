/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
