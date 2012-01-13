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


package org.broad.igv.feature;


import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.WindowFunction;

//~--- JDK imports ------------------------------------------------------------


/**
 * A sub region of a feature.  For example,  a Gene exon
 *
 * @author jrobinso
 */
public class Exon extends AbstractFeature {

    /**
     * The index of the exon relative to the start codon.  The exon with the start
     * codon is number "1".
     */
    private int number;
    private int readingFrame = -1;

    /**
     * Coding start position.  This is the leftmost position of the coding region, not neccessarily the 5'utr end
     */
    private int codingStart;
    private int codingEnd;
    private AminoAcidSequence aminoAcidSequence;
    boolean utr = false;

    // The position of the first base of this exon relative to the start of the mRNA.  This will correspond
    // to either the beginning or end of the exon, depending on the strand
    private int mrnaBase = -1;


    public void setMrnaBase(int base) {
        this.mrnaBase = base;
    }

    public int getAminoAcidNumber(int genomeCoordinate) {
        if (mrnaBase < 0) {
            return -1;
        }
        if (genomeCoordinate < getStart() || genomeCoordinate > getEnd()) {
            throw new IndexOutOfBoundsException();
        }
        if (getStrand() == Strand.POSITIVE) {
            int mrnaCoord = mrnaBase + (genomeCoordinate - codingStart) - 1;
            return mrnaCoord < 0 ? -1 : mrnaCoord / 3 + 1;

        } else if (getStrand() == Strand.NEGATIVE) {
            int mrnaCoord = mrnaBase + (codingEnd - genomeCoordinate);
            return mrnaCoord < 0 ? -1 : mrnaCoord / 3 + 1;

        } else {
            return 0;
        }
    }


    public Exon(String chr, int start, int end, Strand strand) {
        super(chr, start, end, strand);

        // By default the entre exon is a coding region
        this.codingStart = start;
        this.codingEnd = end;
    }

    /**
     * Flag indicating that the entire exon is the UTR.
     *
     * @param utr
     */
    public void setUTR(boolean utr) {
        this.utr = utr;
        if (utr) {
            if (getStrand() == Strand.POSITIVE) {
                codingStart = codingEnd = getEnd();
            } else {
                codingStart = codingEnd = getStart();
            }
        }
    }

    public boolean isUTR() {
        return utr;
    }


    public boolean isUTR(int position) {
        return utr || (position < codingStart || position > codingEnd);
    }

    public void setCodingStart(int codingStart) {
        this.codingStart = Math.max(getStart(), codingStart);
    }


    public void setCodingEnd(int codingEnd) {
        this.codingEnd = Math.min(getEnd(), codingEnd);
    }


    public void setReadingFrame(int offset) {
        this.readingFrame = offset;
    }


    public void setPhase(int phase) {
        if (getStrand() == Strand.POSITIVE) {
            readingFrame = phase;
        } else if (getStrand() == Strand.NEGATIVE) {
            int modLen = (getCodingLength() - phase) % 3;
            readingFrame = modLen;
        }
    }


    public int getCdStart() {
        return codingStart;
    }


    public int getCdEnd() {
        return this.codingEnd;
    }


    public int getCodingLength() {
        return utr ? 0 : Math.max(0, codingEnd - codingStart);
    }


    public int getReadingShift() {
        return readingFrame;
    }

    public AminoAcidSequence getAminoAcidSequence(Genome genome) {
        if (aminoAcidSequence == null) {
            computeAminoAcidSequence(genome);
        }
        return aminoAcidSequence;
    }

    private void computeAminoAcidSequence(Genome genome) {
        if (utr) {
            return;
        }
        int start = getStart();
        int end = getEnd();
        String chr = getChr();
        if (readingFrame >= 0) {
            int readStart = (codingStart > start) ? codingStart : start + readingFrame;
            int readEnd = Math.min(end, codingEnd);
            if (readEnd > readStart + 3) {
                byte[] seqBytes = genome.getSequence(chr, readStart, readEnd);
                aminoAcidSequence = AminoAcidManager.getAminoAcidSequence(seqBytes, readStart, getStrand());
            }
        }
    }


//    public AminoAcid getAminoAcid(Genome genome, int position) {
//        if (aminoAcidSequence == null) {
//            computeAminoAcidSequence(genome);
//        }
//        if (aminoAcidSequence == null) {
//            return null;
//        }
//
//        int delta = position - aminoAcidSequence.getStartPosition();
//        if (delta < 0) {
//            System.out.println("Warning: negative amino acid offset");
//        }
//        int aaIndex = delta / 3;
//        return aaIndex >= aminoAcidSequence.getSequence().size() ? null : aminoAcidSequence.getSequence().get(aaIndex);
//
//    }

    public LocusScore copy() {
        Exon copy = new Exon(getChr(), getStart(), getEnd(), getStrand());
        copy.aminoAcidSequence = this.aminoAcidSequence;
        copy.codingEnd = this.codingEnd;
        copy.codingStart = this.codingStart;
        return copy;
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        String msg = number > 0 ? "Exon number: " + number : "";
        int aaNumber = this.getAminoAcidNumber((int) position);
        if (aaNumber > 0) {
            msg += "<br>Amino acid number: " + aaNumber;
        }
        msg += "<br>" + getLocusString();
        if (description != null) msg += description;
        return msg;
    }

    public void setNumber(int number) {
        this.number = number;
    }

    public String getURL() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

}
