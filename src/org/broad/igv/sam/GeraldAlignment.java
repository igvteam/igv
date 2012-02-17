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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public class GeraldAlignment extends AbstractAlignment implements Alignment {

    private boolean passedFilter;
    private String readSequence = null;
    private int start;
    private int end;

    public GeraldAlignment(String name) {
        this.readName = name;
    }

    public void setReads(String chr, int start, byte[] reads, byte[] qualities) {
        this.chr = chr;
        this.insertions = new AlignmentBlock[0];
        this.alignmentBlocks = new AlignmentBlock[1];
        this.alignmentBlocks[0] = new AlignmentBlock(start, reads, qualities, this);
        this.start = start;
        this.end = start + reads.length;
    }

    public String getReadSequence() {
        if (readSequence == null) {
            readSequence = new String(this.alignmentBlocks[0].getBases());
        }
        return readSequence;
    }

    @Override
    public String getValueString(double position, WindowFunction windowFunction) {
        String str = super.getValueString(position, windowFunction);
        if (!passedFilter) {
            str += "<br>--------------<br>" + "FAILED QUALITY FILTERING";
        }
        return str;
    }

    public String getCigarString() {
        return "*";
    }

    public AlignmentBlock[] getInsertions() {
        return insertions;
    }

    public boolean isProperPair() {
        return isPaired();
    }

    public int getAlignmentStart() {
        return alignmentBlocks[0].getStart();
    }

    public boolean isDuplicate() {
        return false;
    }

    public boolean isMapped() {
        return getChr() != null;
    }

    public boolean isPaired() {
        return this.getMate() != null;
    }

    /**
     * @return the passedFilter
     */
    public boolean isPassedFilter() {
        return passedFilter;
    }

    /**
     * @param passedFilter the passedFilter to set
     */
    public void setPassedFilter(boolean passedFilter) {
        this.passedFilter = passedFilter;
    }

    public int getAlignmentEnd() {
        return end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getSample() {
        return null;
    }

	public boolean isFirstOfPair() {
		return false;
	}

	public boolean isSecondOfPair() {
		return false;
	}

    public Strand getFirstOfPairStrand() {
        return Strand.NONE;
    }

    public Strand getSecondOfPairStrand() {
       return Strand.NONE;
    }

}
