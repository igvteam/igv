package org.igv.sam;

import org.igv.feature.LocusScore;
import org.igv.feature.Strand;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Jan 26, 2011
 */
public class PairedAlignment implements Alignment {

    String chr;
    int start;
    int end;
    Alignment firstAlignment;
    Alignment secondAlignment;



    public PairedAlignment(Alignment firstAlignment) {
        this.firstAlignment = firstAlignment;
        this.start = firstAlignment.getStart();
        this.end = firstAlignment.getEnd();
        this.chr = firstAlignment.getChr();
    }

    public void setSecondAlignment(Alignment alignment) {

        secondAlignment = alignment;
        end = secondAlignment.getEnd();
        // TODO -- check the chrs are equal,  error otherwise

    }

    public String getReadName() {
        return firstAlignment == null ? "" : firstAlignment.getReadName();
    }

    public String getChromosome() {
        return chr;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getChr() {
        return chr;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getContig() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {

        return end;
    }

    public int getAlignmentStart() {
        return start;
    }


    public boolean isMapped() {
        return true;
    }


    public AlignmentBlock[] getAlignmentBlocks() {
        return new AlignmentBlockImpl[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    AlignmentBlock[] insertions;

    public AlignmentBlock[] getInsertions() {
        if (insertions == null) {
            AlignmentBlock[] block1 = firstAlignment.getInsertions();
            if (secondAlignment == null) {
                insertions = block1;
            } else {
                AlignmentBlock[] block2 = secondAlignment.getInsertions();
                insertions = new AlignmentBlockImpl[block1.length + block2.length];
                System.arraycopy(block1, 0, insertions, 0, block1.length);
                System.arraycopy(block2, 0, insertions, block1.length, block2.length);
            }
        }
        return insertions;
    }

    @Override
    public List<Gap> getGaps() {
        return null;
    }

    public int getInferredInsertSize() {
        return Math.abs(firstAlignment.getInferredInsertSize());
    }


    public int getMappingQuality() {
        if (secondAlignment == null) {
            return firstAlignment.getMappingQuality();
        } else {
            return Math.max(firstAlignment.getMappingQuality(), secondAlignment.getMappingQuality());
        }
    }


    public boolean isDuplicate() {
        return firstAlignment.isDuplicate() &&
                (secondAlignment == null || secondAlignment.isDuplicate());
    }


    public boolean isProperPair() {
        return true;  //
    }

    public int getAlignmentEnd() {
        return end;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public byte getBase(double position) {
        if (firstAlignment.contains(position)) {
            return firstAlignment.getBase(position);
        } else if (secondAlignment != null && secondAlignment.contains(position)) {
            return secondAlignment.getBase(position);
        }
        return 0;
    }


    public byte getPhred(double position) {
        if (firstAlignment.contains(position)) {
            return firstAlignment.getPhred(position);
        } else if (secondAlignment != null && secondAlignment.contains(position)) {
            return secondAlignment.getPhred(position);
        }
        return 0;
    }

    /**
     * Return a string to be used for popup text.   The WindowFunction is passed
     * in so it can be used t annotate the value.  The LocusScore object itself
     * does not "know" from what window function it was derived
     *
     * @param mouseX
     * @param renderOptions
     * @return
     */
    public String getAlignmentValueString(double position, int mouseX, AlignmentTrack.RenderOptions renderOptions) {
        StringBuffer buf = new StringBuffer();
        if (secondAlignment != null) {
            buf.append("<table><tr><td valign=\"top\">");
        }
        buf.append("<b>Left alignment</b><br/>");
        buf.append(firstAlignment.getAlignmentValueString(position, mouseX, renderOptions));
        if (secondAlignment != null) {
            buf.append("</td><td valign=\"top\">");
            buf.append("<b>Right alignment</b><br/>");
            buf.append(secondAlignment.getAlignmentValueString(position, mouseX, renderOptions));
            buf.append("</td></tr></table>");
        }
        return buf.toString();
    }

    public boolean contains(double location) {
        return location >= start && location <= end;
    }

    public String getReadSequence() {
        return firstAlignment == null ? "" : firstAlignment.getReadSequence();
    }

    @Override
    public boolean isPaired() {
        return false;  //Counter intuitive, but the pair does not have a mate
    }

    public boolean isNegativeStrand() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setMateSequence(String sequence) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getPairOrientation() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isSmallInsert() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isVendorFailedRead() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Color getYcColor() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getLibrary() {
        return firstAlignment.getLibrary();
    }

    public void setStart(int start) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setEnd(int end) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public float getScore() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public LocusScore copy() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public Alignment getFirstAlignment() {
        return firstAlignment;
    }

    public Alignment getSecondAlignment() {
        return secondAlignment;
    }

    public boolean isFirstOfPair() {
        return false;
    }

    public boolean isSecondOfPair() {
        return false;
    }

    public Strand getFirstOfPairStrand() {
        return firstAlignment.getFirstOfPairStrand();
    }

    public Strand getSecondOfPairStrand() {
        return firstAlignment.getSecondOfPairStrand();
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public void finish() {
        firstAlignment.finish();
        if(secondAlignment != null) {
            secondAlignment.finish();
        }
    }

    @Override
    public boolean isPrimary() {
        return firstAlignment.isPrimary() && (secondAlignment == null || secondAlignment.isPrimary());
    }

    @Override
    public boolean isSupplementary() {
        return firstAlignment.isSupplementary() && (secondAlignment == null || secondAlignment.isSupplementary());
    }

    /**
     *
     * @return the individual Alignment at location or null if neither contains it, if they overlap first is preferred
     */
    @Override
    public Alignment getSpecificAlignment(final double location) {
        Alignment first = getFirstAlignment();
        Alignment second = getSecondAlignment();
        if (first.contains(location)) {
            return first;
        } else if (second.contains(location)) {
            return second;
        } else {
            return null;
        }
    }
}
