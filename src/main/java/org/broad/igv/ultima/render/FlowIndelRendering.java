package org.broad.igv.ultima.render;

import htsjdk.samtools.SAMRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.*;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ultima.FlowUtil;

import java.awt.*;

public class FlowIndelRendering {

    private static ColorMap indelColorMap = ColorMap.getJet(42);

    public boolean handlesAlignment(final Alignment alignment) {

        if ( !(alignment instanceof SAMAlignment) )
            return false;
        final SAMAlignment samAlignment = (SAMAlignment)alignment;

        // must be a flow
        return FlowUtil.isFlow(samAlignment.getRecord());
    }

    public void renderSmallInsertion(Alignment alignment,
                                     AlignmentBlock aBlock,
                                     RenderContext context,
                                     int h, int x, int y,
                                     AlignmentTrack.RenderOptions renderOptions) {

        int pxWing = (h > 10 ? 2 : (h > 5) ? 1 : 0);
        int hairline = 2;
        if ( renderOptions.isInsertQualColoring() ) {
            pxWing = Math.min(pxWing, Math.max(1, (int) (1 / context.getScale())));
            hairline = Math.min(hairline, pxWing);
        }
        Graphics2D g = context.getGraphics();
        g.fillRect(x, y, hairline, h);
        g.fillRect(x - pxWing, y, hairline + 2 * pxWing, hairline);
        g.fillRect(x - pxWing, y + h - hairline, hairline + 2 * pxWing, hairline);

        if ( renderOptions.isInsertQualColoring() && FlowUtil.isFlow(alignment) ) {
            // Ultima: single (==1) INSERT case
            // map qual into a sort of a linear scale and add indicator
            double p = qualsAsProbInsertTP((SAMAlignment) alignment, aBlock);
            if ( p != 0 ) {
                double q = -10 * Math.log10(p);
                Color qColor = new Color(indelColorMap.getColor((int) q));
                g.setColor(qColor);
                g.fillRect(x - pxWing, (int) (y + (h - hairline) * (q / 42)) - 1, hairline + 2 * pxWing, hairline * 2);
            }
        }
    }

    public void renderSmallInsertionWings(Alignment alignment,
                                     AlignmentBlock insertionBlock,
                                     RenderContext context,
                                     int pxH, int pxTop, int pxRight, int pxLeft,
                                     AlignmentTrack.RenderOptions renderOptions) {

        int pxWing = (pxH > 10 ? 2 : 1);  // width of the cursor "wing"
        Graphics2D g = context.getGraphics();

        // adjust wing and hairline
        int hairline = 2;
        double locScale = context.getScale();
        if ( renderOptions.isInsertQualColoring() ) {
            pxWing = Math.min(pxWing, Math.max(1, (int) (1 / locScale)));
            hairline = Math.min(hairline, pxWing);
        }

        g.fillRect(pxLeft - pxWing, pxTop, pxRight - pxLeft + hairline * pxWing, hairline);
        g.fillRect(pxLeft - pxWing, pxTop + pxH - hairline, pxRight - pxLeft + hairline * pxWing, hairline);
        if ( renderOptions.isInsertQualColoring() && FlowUtil.isFlow(alignment) ) {
            // Ultima: large (>1) INSERT case
            // map qual into a sort of a linear scale and add indicator
            double p = qualsAsProbInsertTP((SAMAlignment)alignment, insertionBlock);
            if ( p != 0 ) {
                double q = -10 * Math.log10(p);
                Color qColor = new Color(indelColorMap.getColor((int) q));
                g.setColor(qColor);
                g.fillRect(pxLeft - pxWing, (int) (pxTop + (pxH - hairline) * (q / 42)), pxRight - pxLeft + hairline * pxWing, hairline);
            }
        }

    }

    public void renderDeletionGap(Alignment alignment,
                                          Gap gap,
                                          int y, int h, int x, int w,
                                          RenderContext context,
                                          AlignmentTrack.RenderOptions renderOptions) {

        // collect quals (experimental)
        Color[]       markerColor = new Color[2];
        double[]      markerQ = new double[2];
        if ( renderOptions.isInsertQualColoring() && FlowUtil.isFlow(alignment) ) {

            // locate block who's end is the same as the start of the gap
            boolean blockFound = false;
            int blockIndex = 0;
            final AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
            for ( ; blockIndex < blocks.length ; blockIndex++ ) {
                if ( blocks[blockIndex].getEnd() == gap.getStart() ) {
                    blockFound = true;
                    break;
                }
                else if ( blocks[blockIndex].getEnd() > gap.getStart() )
                    break;
            }
            if ( blockFound && (blockIndex < blocks.length - 1) ) {
                AlignmentBlock abPrev = alignment.getAlignmentBlocks()[blockIndex];
                AlignmentBlock abNext = alignment.getAlignmentBlocks()[blockIndex + 1];
                if (abPrev.getQualities().length > 0 && abNext.getQualities().length > 0) {

                    // calc based on reference
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
                    char        gapBase0 = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart()));
                    char        gapBase1 = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases() - 1));
                    char        alignBase0 = Character.toUpperCase((char)abPrev.getBases().getByte(abPrev.getBases().length - 1));
                    char        alignBase1 = Character.toUpperCase((char)abNext.getBases().getByte(0));
                    char        gapBase0p = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() - 1));
                    char        gapBase1n = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases()));
                    boolean     snp0 = gapBase0p != alignBase0;
                    boolean     snp1 = gapBase1n != alignBase1;
                    byte[]      quals01 = new byte[1];
                    if ( !snp0 ) {
                        quals01[0] = abPrev.getQualities().getByte(abPrev.getQualities().length - 1);
                        double  p;
                        int         delLength = 1;
                        while ( (delLength + 1) < gap.getnBases() &&
                                (gapBase0 == Character.toUpperCase(genome.getReference(alignment.getChr(), gap.getStart() + delLength))) )
                            delLength++;

                        p = qualsAsProbDeleteTP(((SAMAlignment) alignment), abPrev, delLength, false, gap);
                        if ( p != 0 ) {
                            markerQ[0] = -10 * Math.log10(p);
                            markerColor[0] = new Color(indelColorMap.getColor((int) markerQ[0]));
                        }
                    }
                    if ( !snp1 ) {
                        quals01[0] = abNext.getQualities().getByte(0);
                        double  p;
                        int         delLength = 1;
                        while ( (delLength + 1) < gap.getnBases() &&
                                (gapBase1 == Character.toUpperCase(genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases() - delLength))) )
                            delLength++;

                        p = qualsAsProbDeleteTP((SAMAlignment) alignment, abNext, delLength, true, gap);
                        if ( p != 0 ) {
                            markerQ[1] = -10 * Math.log10(p);
                            markerColor[1] = new Color(indelColorMap.getColor((int) markerQ[1]));
                        }
                    }
                }
            }

            // draw delete markers
            Graphics2D g = context.getGraphics();
            if ( markerColor[0] != null || markerColor[1] != null ) {

                int hairline = Math.min(2, (int) (1 / context.getScale()));

                Color c = g.getColor();
                if ((markerQ[0] == markerQ[1]) || (gap.getnBases() == 1)) {

                    // draw a full line, average as needed
                    double q1 = markerQ[0];
                    Color c1 = markerColor[0];
                    if (c1 == null) {
                        q1 = markerQ[1];
                        c1 = markerColor[1];
                    } else if (markerColor[1] != null && markerQ[1] != q1) {
                        q1 = (q1 + markerQ[1]) / 2;
                        c1 = new Color(indelColorMap.getColor((int) q1));
                    }
                    g.setColor(c1);
                    g.fillRect(x + (int) (w * q1 / 42) - hairline, y, hairline * 2, h);
                } else {
                    if (markerColor[0] != null) {
                        g.setColor(markerColor[0]);
                        g.fillRect(x + (int) (w * markerQ[0] / 42) - hairline, y, hairline * 2, h * 3 / 4);
                    }
                    if (markerColor[1] != null) {
                        g.setColor(markerColor[1]);
                        g.fillRect(x + (int) (w * markerQ[1] / 42) - hairline, y + h / 4, hairline * 2, h * 3 / 4);
                    }
                }

                g.setColor(c);
            }
        }
    }

    static private double qualsAsProb(ByteSubarray quals) {

        // calc prob
        double              probSum = 0.0;
        double              probCount = 0;
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final byte q = quals.getByte(i);
            if (q != 255) {
                probSum += Math.pow(10.0, -q / 10.0);
                probCount++;

            }
        }
        if ( probCount != 0 )
            return probSum / probCount;
        else
            return -1;
    }


    public double qualsAsProbInsertTP(SAMAlignment samAlignment, AlignmentBlock block) {
        return qualsAsProbInsertTP(samAlignment, block, 0, block.getLength());
    }

    private class HMer {
        int     start;
        int     end;
    }

    public double qualsAsProbInsertTP(SAMAlignment samAlignment, AlignmentBlock block, int fragOfs, int fragLength) {

        // for now, handle only hmer inserts
        if ( !blockIsHmer(block, fragOfs, fragLength) ) {

            // by definition, blocks failing this test can only be of multiple bases
            // try breaking them into two block fragments, one from front and one from back
            byte[]  bases = block.getBases().getBytes();
            int     f1ofs = 0;
            int     f1Length = 0;
            int     base = bases[f1ofs];
            for ( int n = 0 ; n < bases.length ; n++ ) {
                if (bases[n] == base)
                    f1Length++;
                else
                    break;
            }

            int     f2ofs = bases.length;
            int     f2Length = 0;
            base = bases[f2ofs-1];
            for ( int n = bases.length - 1 ; n >= 0 ; n-- ) {
                if (bases[n] == base) {
                    f2Length++;
                    f2ofs--;
                }
                else
                    break;
            }

            double      p1 = qualsAsProbInsertTP(samAlignment, block, f1ofs, f1Length);
            double      p2 = qualsAsProbInsertTP(samAlignment, block, f2ofs, f2Length);

            if ( p1 == 0 )
                return p2;
            else if ( p2 == 0 )
                return p1;
            else
                return (p1 + p2) / 2;
        }

        // access read/record
        SAMRecord record = samAlignment.getRecord();
        if ( record == null )
            return 0;


        // establish the hmer on which this block sits (inside the read...).
        byte        base = block.getBase(fragOfs);
        int         start = block.getIndexOnRead() + fragOfs;

        // short-circuit a simple and common case: insert of 1 with tp=-1 -> quality is already here!
        if ( fragLength == 1 ) {
            byte[]      tp = record.getByteArrayAttribute(FlowUtil.ATTR_TP);
            if ( tp[start] == -1 ) {
                byte      q = block.getQuality(0);
                return Math.pow(10.0, -q / 10.0);
            }
        }
        HMer        hmer = findHmer(record, start, fragLength, base);

        // find tp value and return it
        return findQualByTPValue(record, hmer, -fragLength);
    }

    private HMer findHmer(SAMRecord record, int start, int length, byte base) {

        HMer        hmer = new HMer();
        hmer.start = start;

        byte[]      bases = record.getReadBases();
        hmer.end = hmer.start + (length - 1);
        while ( hmer.start > 0 && bases[hmer.start - 1] == base )
            hmer.start--;
        while ( (hmer.end + 1) < bases.length  && bases[hmer.end + 1] == base )
            hmer.end++;

        return hmer;
    }

    private double findQualByTPValue(SAMRecord record, HMer hmer, int tpValue) {

        // get quals and tp
        byte[]      tp = record.getByteArrayAttribute(FlowUtil.ATTR_TP);

        // scan for tpValue, extract qual
        for ( int ofs = hmer.start ; ofs <= hmer.end ; ofs++ )
            if ( tp[ofs] == tpValue ) {
                byte[]      quals = record.getBaseQualities();
                return Math.pow(10.0, -quals[ofs] / 10.0);
            }

        // if here, none of the tp values matched
        // check for the special case of a delete which is of an original hmer larger than mc (def:12)
        if ( tpValue > 0 ) {
            int     mc = getMC(record);
            int     hmerSize = hmer.end - hmer.start + 1;
            if ( tpValue + hmerSize > mc )
                return 1.0 - FlowUtil.MIN_PROB_DEFAULT;
        }

        // if here, simply fail
        return 0;
    }

    private int getMC(SAMRecord record) {
        try {
            return Integer.parseInt(record.getReadGroup().getAttribute("mc"));
        } catch (Exception e) {
            return 0;
        }
    }

    private boolean blockIsHmer(AlignmentBlock block, int fragOfs, int fragLength) {
        byte[]      bases = block.getBases().getBytes();

        if ( bases == null || bases.length < (fragOfs + fragLength) )
            return false;
        else if ( fragLength == 1 )
            return true;

        byte base = bases[fragOfs];
        for ( int n = fragOfs + 1 ; n < fragOfs + fragLength ; n++ ) {
            if (bases[n] != base)
                return false;
        }
        return true;
    }

    public double qualsAsProbDeleteTP(SAMAlignment samAlignment, AlignmentBlock block, int delLength, boolean delIsBeforeBlock, Gap gap)
    {
        // access read/record
        SAMRecord   record = samAlignment.getRecord();
        if ( record == null )
            return 0;

        // establish the hmer on which this block sits (inside the read...)
        byte        base = block.getBase(delIsBeforeBlock ? 0 : block.getLength() - 1);
        int         start = block.getIndexOnRead();
        if ( !delIsBeforeBlock )
            start += (block.getLength() - 1);
        HMer        hmer = findHmer(record, start, 0, base);

        // try establising by using t0
        final double t0qual = qualsAsProbDeleteTPByT0(samAlignment, record, block, delLength, delIsBeforeBlock, gap);
        if ( t0qual != 0 ) {
            return t0qual;
        }

        return findQualByTPValue(record, hmer, delLength);
    }

    private double qualsAsProbDeleteTPByT0(SAMAlignment samAlignment, SAMRecord record, AlignmentBlock block, int delLength, boolean delIsBeforeBlock, Gap gap) {

        // consider using t0 only if DEL of one base (additional conditions to follow)
        if ( delLength != 1 ) {
            return 0;
        }

        // get location just before this read
        final int         loc = block.getIndexOnRead() + (delIsBeforeBlock ? -1 : (block.getLength() - 1));
        if ( loc < 0 )
            return 0;

        // get bases around the deletion. They must be different than the deleted base
        final byte[]      bases = record.getReadBases();
        if ( loc >= (bases.length - 1) )
            return 0;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        final byte delBase = (byte) Character.toUpperCase(genome.getReference(samAlignment.getChr(), gap.getStart()));
        if ( (delBase == bases[loc]) ||  (delBase == bases[loc+1]) )
            return 0;

        // extract t0 values for the surrounding bases, convert to prob
        if ( !record.hasAttribute(FlowUtil.TAG_T0) )
            return 0;
        final byte[]      t0 = record.getStringAttribute(FlowUtil.TAG_T0).getBytes();
        final double      p0 = Math.pow(10.0, (t0[loc] - '!') / -10.0);
        final double      p1 = Math.pow(10.0, (t0[loc+1] - '!') / -10.0);

        // return maximal of the two
        return Math.max(p0, p1);
    }
}
