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

package org.broad.igv.feature;


import com.google.common.base.Objects;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.WindowFunction;

import java.lang.reflect.InvocationHandler;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;

//~--- JDK imports ------------------------------------------------------------


/**
 * A sub region of a feature.  For example,  a Gene exon
 *
 * @author jrobinso
 */
public class Exon extends AbstractFeature implements IExon {

    /**
     * Index relative to the 5' end.
     */
    private int number;

    /**
     * Coding start position.  This is the leftmost position of the coding region, not necessarily the 5'utr end
     */
    private int codingStart;
    private int codingEnd;

    private AminoAcidSequence aminoAcidSequence;

    //Raw bytes representing nucleotides
    //Stored separately from the aminoAcidSequence because the latter changes
    //when we change translation tables
    private byte[] seqBytes;

    private boolean noncoding = false;

    // The position of the first base of this exon relative to the start of the mRNA.  This will correspond
    // to either the beginning or end of the exon, depending on the strand
    private int mrnaBase = -1;

    public void setMrnaBase(int base) {
        this.mrnaBase = base;
    }

    /**
     * Get amino acid number based on genomic coordinate.
     * Genome coordinate MUST be 0-based
     *
     * @param genomeCoordinate
     * @return
     */
    public int getAminoAcidNumber(int genomeCoordinate) {
        if (mrnaBase < 0) {
            return -1;
        }
        if (genomeCoordinate < getStart() || genomeCoordinate > getEnd()) {
            throw new IndexOutOfBoundsException();
        }
        if (getStrand() == Strand.POSITIVE) {
            int mrnaCoord = mrnaBase + (genomeCoordinate - codingStart);
            return mrnaCoord < 0 ? -1 : mrnaCoord / 3 + 1;

        } else if (getStrand() == Strand.NEGATIVE) {
            //Since codingEnd is exclusive-end, we subtract 1 from it
            //We want mrnaCoord = 0 when genomeCoordinate == codingEnd - 1
            int mrnaCoord = mrnaBase + (codingEnd - 1 - genomeCoordinate);
            return mrnaCoord < 0 ? -1 : mrnaCoord / 3 + 1;

        } else {
            return 0;
        }
    }


    public Exon(String chr, int start, int end, Strand strand) {
        super(chr, start, end, strand);

        // By default the entire exon is a coding region
        this.codingStart = start;
        this.codingEnd = end;
    }


    public Exon(Exon bf) {
        this.start = bf.getStart();
        this.end = bf.getEnd();
        this.strand = bf.getStrand();
        this.codingStart = bf.getCdStart();
        this.codingEnd = bf.getCdEnd();
        this.chromosome = bf.getChr();
        this.type = bf.getType();
        this.color = bf.getColor();
        this.description = bf.getDescription();
        this.attributes = bf.getAttributes();
        this.name = bf.getName();
        this.readingFrame = bf.getReadingFrame();

        this.noncoding = (type != null && SequenceOntology.utrTypes.contains(type));
    }


    public Exon(BasicFeature bf) {
        this.start = bf.getStart();
        this.end = bf.getEnd();
        this.strand = bf.getStrand();
        this.codingStart = bf.getThickStart();
        this.codingEnd = bf.getThickEnd();
        this.chromosome = bf.getChr();
        this.type = bf.getType();
        this.color = bf.getColor();
        this.description = bf.getDescription();
        this.attributes = bf.getAttributes();
        this.name = bf.getName();
        this.readingFrame = bf.getReadingFrame();

        this.noncoding = (type != null && SequenceOntology.utrTypes.contains(type));
    }

    /**
     * Flag indicating that the entire exon is non-coding.
     *
     * @param bool
     */
    public void setNonCoding(boolean bool) {
        this.noncoding = bool;
        if (bool) {
            if (getStrand() == Strand.POSITIVE) {
                codingStart = codingEnd = getEnd();
            } else {
                codingStart = codingEnd = getStart();
            }
        }
    }

    public boolean isNonCoding() {
        return noncoding;
    }

    public void setCodingStart(int codingStart) {
        this.codingStart = Math.max(getStart(), codingStart);
    }


    public void setCodingEnd(int codingEnd) {
        this.codingEnd = Math.min(getEnd(), codingEnd);
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
        return noncoding ? 0 : Math.max(0, codingEnd - codingStart);
    }

    public AminoAcidSequence getAminoAcidSequence(Genome genome, Exon prevExon, Exon nextExon) {
        if (aminoAcidSequence == null ||
                //If the stored sequence was computed with a different codon table, we reset
                !(Objects.equal(aminoAcidSequence.getCodonTableKey(), AminoAcidManager.getInstance().getCodonTable().getKey()))) {
            computeAminoAcidSequence(genome, prevExon, nextExon);
        }
        return aminoAcidSequence;
    }

    private void computeAminoAcidSequence(Genome genome, Exon prevExon, Exon nextExon) {
        if (noncoding) {
            return;
        }
        int start = getStart();
        int end = getEnd();
        String chr = getChr();
        if (readingFrame >= 0) {
            int readStart = (codingStart > start) ? codingStart : start;
            int readEnd = Math.min(end, codingEnd);
            if (readEnd > readStart + 3) {
                if (seqBytes == null) {
                    seqBytes = genome.getSequence(chr, readStart, readEnd);
                }
                if (seqBytes != null) {

                    // Grab nucleotides from previous exon if needed to complete first codon
                    if (readingFrame > 0 && prevExon != null) {
                        int diff = 3 - readingFrame;
                        byte[] d = genome.getSequence(chr, prevExon.getCdEnd() - diff, prevExon.getCdEnd());
                        byte[] tmp = new byte[d.length + seqBytes.length];
                        System.arraycopy(d, 0, tmp, 0, diff);
                        System.arraycopy(seqBytes, 0, tmp, diff, seqBytes.length);
                        seqBytes = tmp;
                        readStart -= diff;
                    }

                    // Grab nucleotides from next exon if needed for last codon
                    int diff = 3 - ((readEnd - (codingStart + readingFrame)) % 3);
                    if (diff > 0 && diff < 3 && nextExon != null) {
                        byte[] d = genome.getSequence(chr, nextExon.getCdStart(), nextExon.getCdStart() + diff);
                        byte[] tmp = new byte[d.length + seqBytes.length];
                        System.arraycopy(seqBytes, 0, tmp, 0, seqBytes.length);
                        System.arraycopy(d, 0, tmp, seqBytes.length, d.length);
                        seqBytes = tmp;
                    }

                    aminoAcidSequence = AminoAcidManager.getInstance().getAminoAcidSequence(getStrand(), readStart, seqBytes);
                }
            }
        }
    }

    public Exon copy() {
        Exon copy = new Exon(getChr(), getStart(), getEnd(), getStrand());
        copy.seqBytes = this.seqBytes;
        copy.aminoAcidSequence = this.aminoAcidSequence;
        copy.codingEnd = this.codingEnd;
        copy.codingStart = this.codingStart;
        copy.name = this.name;
        copy.noncoding = this.noncoding;
        copy.mrnaBase = this.mrnaBase;
        return copy;
    }

    public String getValueString(double position, WindowFunction windowFunction) {

        StringBuffer buffer = new StringBuffer();

        if (number > 0) buffer.append("Exon number: " + number + "<br>");
        int aaNumber = this.getAminoAcidNumber((int) position);
        if (aaNumber > 0) {
            buffer.append("Amino acid codingNumber: " + aaNumber + "<br>");
        }
        buffer.append(getLocusString());
        if (description != null) buffer.append("<br>" + description);
        if (attributes != null) {
            buffer.append(getAttributeString());
        }

        return buffer.toString();
    }

    public void setNumber(int number) {
        this.number = number;
    }

    public String getURL() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public static IExon getExonProxy(IExon exon) {
        InvocationHandler handler = new ExonLocHandler(exon);
        IExon eProx = (IExon) Proxy.newProxyInstance(IExon.class.getClassLoader(),
                new Class[]{IExon.class},
                handler);
        return eProx;
    }

    private static class ExonLocHandler implements InvocationHandler {

        private IExon parent;
        private int hashCode = 0;

        public ExonLocHandler(IExon parent) {
            this.parent = parent;
        }

        private boolean equals(IExon parent, Object inother) {
            if (inother == null || !(inother instanceof IExon)) {
                return false;
            }
            IExon other = (IExon) inother;
            boolean eq = parent.getChr().equals(other.getChr());
            eq &= parent.getStart() == other.getStart();
            eq &= parent.getEnd() == other.getEnd();
            eq &= parent.getCdStart() == other.getCdStart();
            eq &= parent.getCdEnd() == other.getCdEnd();
            eq &= parent.getStrand() == other.getStrand();
            return eq;
        }

        private int hashCode(IExon parent) {
            if (hashCode != 0) {
                return hashCode;
            }

            String conc = parent.getChr() + parent.getStrand().toString() + parent.getStart();
            conc += parent.getEnd();
            conc += parent.getCdStart();
            conc += parent.getCdEnd();
            int hc = conc.hashCode();

            if (hc == 0) {
                hc = 1;
            }
            hashCode = hc;
            return hc;
        }

        @Override
        public Object invoke(Object proxy, Method method, Object[] args) throws Throwable {
            if (method.getName().equals("hashCode")) {
                return hashCode(parent);
            } else if (method.getName().equals("equals")) {
                return equals(parent, args[0]);
            } else {
                return method.invoke(parent, args);
            }
        }
    }

}
