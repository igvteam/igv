/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
        this.alignmentBlocks[0] = new AlignmentBlock(getChr(), start, reads, qualities);
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
