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

import java.awt.*;
import java.text.DecimalFormat;

/**
 * @author jrobinso
 */
public class GeraldAlignment implements Alignment {

    private boolean passedFilter;
    private String readSequence = null;
    private int start;
    private int end;
    String chr;
    int inferredInsertSize;
    int mappingQuality = 255;  // 255 by default
    ReadMate mate;
    String readName;
    AlignmentBlock[] alignmentBlocks;
    AlignmentBlock[] insertions;
    char[] gapTypes;
    private boolean negativeStrand;
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();

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

    @Override
    public String getReadName() {
        return null;
    }

    public String getReadSequence() {
        if (readSequence == null) {
            readSequence = new String(this.alignmentBlocks[0].getBases());
        }
        return readSequence;
    }

    @Override
    public String getChromosome() {
        return null;
    }

    @Override
    public String getChr() {
        return null;
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = null;

        // First check insertions.  Position is zero based, block coords 1 based
        if (this.insertions != null) {
            for (AlignmentBlock block : this.insertions) {
                double insertionLeft = block.getStart() - .25;
                double insertionRight = block.getStart() + .25;
                if (position > insertionLeft && position < insertionRight) {
                    if (block.hasFlowSignals()) {
                        int offset;
                        buf = new StringBuffer();
                        buf.append("Insertion: " + new String(block.getBases()) + "<br>");
                        buf.append("Base phred quality = ");
                        for (offset = 0; offset < block.getLength(); offset++) {
                            byte quality = block.getQuality(offset);
                            if (0 < offset) {
                                buf.append(",");
                            }
                            buf.append(quality);
                        }
                        buf.append("<br>");
                        for (offset = 0; offset < block.getLength(); offset++) {
                            byte base = block.getBase(offset);
                            buf.append((char) base + ": ");
                        }
                        buf.append("----------------------"); // NB: no <br> required
                        return buf.toString();
                    } else {
                        return "Insertion: " + new String(block.getBases());
                    }
                }
            }
        }

        buf = new StringBuffer();

        String sample = getSample();
        if (sample != null) {
            buf.append("Sample = " + sample + "<br>");
        }
        String readGroup = getReadGroup();
        if (sample != null) {
            buf.append("Read group = " + readGroup + "<br>");
        }
        buf.append("----------------------" + "<br>");

        int basePosition = (int) position;
        buf.append("Read name = " + getReadName() + "<br>");
        buf.append("Location = " + getChr() + ":" + DECIMAL_FORMAT.format(1 + (long) position) + "<br>");
        buf.append("Alignment start = " + DECIMAL_FORMAT.format(getAlignmentStart() + 1) + " (" + (isNegativeStrand() ? "-" : "+") + ")<br>");
        buf.append("Cigar = " + getCigarString() + "<br>");
        buf.append("Mapped = " + (isMapped() ? "yes" : "no") + "<br>");
        buf.append("Mapping quality = " + getMappingQuality() + "<br>");
        buf.append("----------------------" + "<br>");

        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte base = block.getBase(offset);
                byte quality = block.getQuality(offset);
                buf.append("Base = " + (char) base + "<br>");
                buf.append("Base phred quality = " + quality + "<br>");
                if (block.hasCounts()) {
                    buf.append("Count = " + block.getCount(offset) + "<br>");
                }
            }
        }

        if (this.isPaired()) {
            buf.append("----------------------" + "<br>");
            buf.append("Mate start = " + getMate().positionString() + "<br>");
            buf.append("Mate is mapped = " + (getMate().isMapped() ? "yes" : "no") + "<br>");
            //buf.append("Pair is proper = " + (getProperPairFlag() ? "yes" : "no") + "<br>");
            if (getChr().equals(getMate().getChr())) {
                buf.append("Insert size = " + getInferredInsertSize() + "<br>");
            }
            if (getPairOrientation().length() > 0) {
                buf.append("Pair orientation = " + getPairOrientation() + "<br>");
            }
        }
        buf.append("----------------------");

        if (!passedFilter) {
            buf.append("<br>--------------<br>" + "FAILED QUALITY FILTERING");
        }

        return buf.toString();
    }


    public String getCigarString() {
        return "*";
    }

    @Override
    public int getInferredInsertSize() {
        return 0;
    }

    @Override
    public int getMappingQuality() {
        return 0;
    }

    @Override
    public ReadMate getMate() {
        return null;
    }

    public AlignmentBlock[] getInsertions() {
        return insertions;
    }

    @Override
    public char[] getGapTypes() {
        return new char[0];
    }

    public boolean isProperPair() {
        return isPaired();
    }

    public int getAlignmentStart() {
        return alignmentBlocks[0].getStart();
    }

    @Override
    public boolean contains(double location) {
        return false;
    }

    @Override
    public AlignmentBlock[] getAlignmentBlocks() {
        return new AlignmentBlock[0];
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

    @Override
    public byte getBase(double position) {
        return 0;
    }

    @Override
    public byte getPhred(double position) {
        return 0;
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




    @Override

    public float getScore() {
        return getMappingQuality();
    }

    public String getSample() {
        return null;
    }

    @Override
    public String getReadGroup() {
        return null;
    }

    @Override
    public Object getAttribute(String key) {
        return null;
    }

    @Override
    public void setMateSequence(String sequence) {

    }

    @Override
    public String getPairOrientation() {
        return "";
    }

    @Override
    public boolean isSmallInsert() {
        int absISize = Math.abs(getInferredInsertSize());
        return absISize > 0 && absISize <= getReadLength();
    }


    public int getReadLength() {
        return getReadSequence().length();
    }

    @Override
    public boolean isVendorFailedRead() {
        return false;
    }

    @Override
    public Color getColor() {
        return null;
    }

    @Override
    public String getLibrary() {
        return null;
    }

    @Override

    public String getClipboardString(double location) {
        return getValueString(location, null);
    }

    @Override

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

    @Override
    public boolean isNegativeStrand() {
        return negativeStrand;
    }

    public void setMappingQuality(int mappingQuality) {
        this.mappingQuality = mappingQuality;
    }

    public void setNegativeStrand(boolean negativeStrand) {
        this.negativeStrand = negativeStrand;
    }

    public void setInferredInsertSize(int inferredInsertSize) {
        this.inferredInsertSize = inferredInsertSize;
    }

    public void setMate(ReadMate mate) {
        this.mate = mate;
    }
}
