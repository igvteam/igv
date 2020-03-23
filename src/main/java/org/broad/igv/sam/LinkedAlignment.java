package org.broad.igv.sam;

import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;
import java.util.*;

/**
 * Class for experimenting with 10X linked reads.
 */

public class LinkedAlignment implements Alignment {


    final String tag;       // Tag used to link, or "READNAME"
    final String name;   // Tag value  (usually barcode or readname)

    String haplotype;
    String sample;
    String readGroup;
    String library;

    String chr;
    Strand strand;
    int alignmentStart;
    int alignmentEnd;

    List<Alignment> alignments;


    Map<String, Object> attributes;

    public LinkedAlignment(String tag, String bc) {
        attributes = new HashMap<>();
        this.tag = tag;
        this.name = bc;
        attributes.put(tag, name);
        alignments = new ArrayList<>();
    }

    public void addAlignment(Alignment alignment) {

        String sample = alignment.getSample();
        String readGroup = alignment.getReadGroup();
        String library = alignment.getLibrary();

        if (alignments.isEmpty()) {
            this.chr = alignment.getChr();
            alignmentStart = alignment.getAlignmentStart();
            alignmentEnd = alignment.getAlignmentEnd();

            Object hp = alignment.getAttribute("HP");
            haplotype = hp == null ? null : hp.toString();

            strand = alignment.getReadStrand();
            this.sample = sample == null ? "" : sample;
            this.readGroup = readGroup == null ? "" : readGroup;
            this.library = library == null ? "" : library;

        } else {

            if (!this.chr.equals(alignment.getChr())) {
                throw new RuntimeException("Mixed chromosome linked alignments not supported");
            }

            alignmentStart = Math.min(alignment.getAlignmentStart(), this.alignmentStart);
            alignmentEnd = Math.max(alignment.getAlignmentEnd(), this.alignmentEnd);

            Object hp = alignment.getAttribute("HP");
            if (hp != null) {
                if (this.haplotype == null) {
                    this.haplotype = hp.toString();
                } else if (!hp.toString().equals(this.haplotype)) {
                    this.haplotype = "MIXED";
                }
            }
            if (this.strand != alignment.getReadStrand()) {
                this.strand = Strand.NONE;   // i.e. mixed
            }
            if (!this.sample.equals(sample)) {
                this.sample += ", " + sample;
            }
            if (!this.readGroup.equals(readGroup)) {
                this.readGroup += ", " + readGroup;
            }
            if (!this.library.equals(library)) {
                this.library += ", " + library;
            }

        }

        alignments.add(alignment);
    }

    public Strand getStrand() {
        return strand;
    }

    /**
     * Return the strand of the linked alignment at the genomic position
     *
     * @param position
     * @return
     */
    public Strand getStrandAtPosition(double position) {
        if (strand == Strand.NONE) {
            for (Alignment a : alignments) {
                if (a.contains(position)) {
                    return a.getReadStrand();
                }
            }
        }
        return strand;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public int getAlignmentStart() {
        return this.alignmentStart;
    }

    @Override
    public int getAlignmentEnd() {
        return this.alignmentEnd;
    }


    @Override
    public int getStart() {
        return this.alignmentStart;
    }

    @Override
    public int getEnd() {
        return this.alignmentEnd;
    }

    @Override
    public boolean contains(double location) {
        return location >= this.alignmentStart && location <= this.alignmentEnd;
    }

    @Override
    public boolean isMapped() {
        return true;
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {

        if (alignments.size() == 1) {
            return alignments.get(0).getValueString(position, mouseX, windowFunction);
        } else {

            // First check to see if we are over an insertion.   Insertions take precedence.
            for (Alignment a : alignments) {
                for (AlignmentBlock block : a.getInsertions()) {
                    if (block.containsPixel(mouseX)) {
                        return a.getValueString(position, mouseX, windowFunction);
                    }
                }
            }

            StringBuffer buffer = new StringBuffer();
            buffer.append("Linking id (" + tag + ") = " + this.name);
            if (this.haplotype != null) buffer.append("<br>Haplotype = " + this.haplotype);
            buffer.append("<br># alignments = " + alignments.size());
            buffer.append("<br>Total span = " + Globals.DECIMAL_FORMAT.format(getAlignmentEnd() - getAlignmentStart()) + "bp");

            // Link by readname == supplementary alignments.  Crude "is not 10x?" test
            if ("READNAME".equals(tag)) {
                buffer.append("<br>Strands = ");
                for (Alignment a : alignments) {
                    buffer.append(a.getReadStrand() == Strand.POSITIVE ? "+" : "-");
                }

                for (Alignment a : alignments) {
                    if (a instanceof SAMAlignment) {
                        buffer.append("<br>");
                        buffer.append(((SAMAlignment) a).getSynopsisString());
                    }
                }
            }

            for (Alignment a : alignments) {
                if (a.contains(position)) {
                    buffer.append("<hr>");
                    buffer.append(a.getValueString(position, mouseX, windowFunction));
                }
            }


            return buffer.toString();
        }
    }

    @Override
    public Object getAttribute(String key) {
        if ("HP".equals(key)) {
            return haplotype;
        } else {
            return attributes.get(key);
        }
    }


    @Override
    public int getMappingQuality() {
        return 30;    // This is used for coloring.  Not sure what to do here
    }

    /////////////////////////////////////////////////////////////

    @Override
    public String getReadName() {
        return "READNAME".equals(tag) ? name : null;
    }

    @Override
    public String getReadSequence() {
        return null;
    }

    @Override
    public AlignmentBlock[] getAlignmentBlocks() {
        return null;
    }

    @Override
    public AlignmentBlock[] getInsertions() {

        int n = 0;
        for (Alignment a : alignments) {
            n += a.getInsertions().length;
        }
        AlignmentBlock[] insertions = new AlignmentBlock[n];

        n = 0;
        for (Alignment a : alignments) {
            AlignmentBlock[] blocks = a.getInsertions();
            System.arraycopy(blocks, 0, insertions, n, blocks.length);
            n += blocks.length;
        }
        return insertions;
    }


    @Override
    public String getCigarString() {
        return null;
    }

    @Override
    public List<Gap> getGaps() {
        return null;
    }

    @Override
    public int getInferredInsertSize() {
        return 0;
    }

    @Override
    public ReadMate getMate() {
        return null;
    }

    @Override
    public Strand getReadStrand() {
        return null;
    }

    @Override
    public boolean isProperPair() {
        return false;
    }


    @Override
    public boolean isPaired() {
        return false;
    }

    @Override
    public boolean isFirstOfPair() {
        return false;
    }

    @Override
    public boolean isSecondOfPair() {
        return false;
    }

    @Override
    public boolean isNegativeStrand() {
        return strand == Strand.NEGATIVE;
    }

    @Override
    public boolean isDuplicate() {
        return false;
    }

    @Override
    public boolean isPrimary() {
        return false;
    }

    @Override
    public boolean isSupplementary() {
        return false;
    }

    @Override
    public byte getBase(double position) {

        byte base = 0;

        for (Alignment al : alignments) {
            if (al.contains(position)) {
                byte b = al.getBase(position);
                if (base == 0) {
                    base = b;
                } else {
                    if (base != b) {
                        base = 0;
                        break;
                    }
                }
            }
        }

        return base;
    }

    @Override
    public byte getPhred(double position) {
        return 0;
    }

    @Override
    public void setMateSequence(String sequence) {

    }

    @Override
    public String getPairOrientation() {
        return null;
    }

    @Override
    public Strand getFirstOfPairStrand() {
        return null;
    }

    @Override
    public Strand getSecondOfPairStrand() {
        return null;
    }

    @Override
    public boolean isVendorFailedRead() {
        return false;
    }

    public Color getYcColor() {
        return null;
    }

    @Override
    public String getSample() {
        return sample;
    }

    @Override
    public String getReadGroup() {

        return null;
    }

    @Override
    public String getLibrary() {
        return null;
    }

    @Override
    public String getClipboardString(double location, int mouseX) {
        return null;
    }

    @Override
    public void finish() {
        alignments.sort(ALIGNMENT_START_COMPARATOR);
    }

    @Override
    public void setStart(int start) {

    }

    @Override
    public void setEnd(int end) {

    }

    @Override
    public float getScore() {
        return 0;
    }

    @Override
    public String getContig() {
        return null;
    }

    static final Comparator<Alignment> ALIGNMENT_START_COMPARATOR = new Comparator<Alignment>() {
        public int compare(Alignment o1, Alignment o2) {
            return o1.getAlignmentStart() - o2.getAlignmentStart();
        }
    };

}
