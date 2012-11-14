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

package org.broad.igv.hic.tools;

/**
 * @author Jim Robinson
 * @date 9/24/11
 */
public class AlignmentPair {

    private boolean strand1;  // true if forward strand
    private boolean strand2;
    private int frag1;
    private int frag2;
    private int chr1;
    private int pos1;
    private int chr2;
    private int pos2;
   /*
    public AlignmentPair(int chr1, int pos1, int chr2, int pos2) {
        this.chr1 = chr1;
        this.pos1 = pos1;
        this.chr2 = chr2;
        this.pos2 = pos2;
    }
     */
    public AlignmentPair(boolean strand1, int chr1, int pos1, int frag1, boolean strand2, int chr2, int pos2, int frag2) {
        this.strand1 = strand1;
        this.chr1 = chr1;
        this.pos1 = pos1;
        this.frag1 = frag1;
        this.strand2 = strand2;
        this.chr2 = chr2;
        this.pos2 = pos2;
        this.frag2 = frag2;
    }


    public int getChr1() {
        return chr1;
    }

    public int getPos1() {
        return pos1;
    }

    public int getChr2() {
        return chr2;
    }

    public int getPos2() {
        return pos2;
    }

    public boolean getStrand1() {
        return strand1;
    }

    public int getStrand1AsInt() {
        return strand1?0:16;
    }

    public int getStrand2AsInt() {
        return strand2?0:16;
    }

    public boolean getStrand2() {
        return strand2;
    }

    public int getFrag1() {
        return frag1;
    }

    public int getFrag2() {
        return frag2;
    }

    public String toString() {
        int str1 = getStrand1AsInt();
        int str2 = getStrand2AsInt();
        return str1 + "\t" +chr1 + "\t" + pos1 + "\t" + frag1 + "\t" + str2 + "\t" + chr2 + "\t" + pos2 +"\t" + frag2;
    }
}
