package org.broad.igv.ultima.render;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.*;
import org.broad.igv.track.RenderContext;

import java.awt.*;

public class FlowIndelRendering {

    // constants
    private static final String TAG_T0 = "t0";
    private static final String ATTR_TP = "tp";
    private static final String RG_ATTR_PL = "PL";
    private static final String RG_ATTR_MC = "mc";
    private static final String RG_ATTR_PL_ULTIMA = "ULTIMA";
    private static final double MIN_POSSIBLE_QUALITY = 0;
    public static final int COLORMAP_SIZE = 44;

    // the color map
    private static ColorMap indelColorMap = ColorMap.getJet(COLORMAP_SIZE);

    // an Hmer
    static class Hmer {
        int start;
        int end;
        int backwardsSize = 0;
        int forwardSize = 0;

        int size() { return end - start + 1; }
    }


    // is this alignment handled by this renderer?
    public boolean handlesAlignment(final Alignment alignment) {

        // we only handle sam alignments
        if ( !(alignment instanceof SAMAlignment) )
            return false;

        return isFlow(((SAMAlignment)alignment).getRecord());
    }

    public static boolean isFlow(SAMRecord record) {

        final SAMReadGroupRecord readGroup = record.getReadGroup();
        if ( readGroup == null )
            return false;
        if ( !RG_ATTR_PL_ULTIMA.equals(readGroup.getAttribute(RG_ATTR_PL))
                &&  (readGroup.getAttribute(RG_ATTR_MC) == null) )
            return false;
        if ( !record.hasAttribute(ATTR_TP)  )
            return false;

        return true;
    }

    // is this block (INS) handled (accepted)
    public boolean handlesBlock(AlignmentBlock block) {
        return true;
    }

    // is this gap (DEL) handled (accepted)
    public boolean handlesGap(Gap gap) {
        return true;
    }

    public void renderSmallInsertion(Alignment alignment,
                                     AlignmentBlock block,
                                     RenderContext context,
                                     int h, int x, int y,
                                     AlignmentTrack.RenderOptions renderOptions) {

        int pxWing = (h > 10 ? 2 : (h > 5) ? 1 : 0);
        int hairline = 2;
        if ( renderOptions.isIndelQualColoring() ) {
            pxWing = Math.min(pxWing, Math.max(1, (int) (1 / context.getScale())));
            hairline = Math.min(hairline, pxWing);
        }
        Graphics2D g = context.getGraphics();
        g.fillRect(x, y, hairline, h);
        g.fillRect(x - pxWing, y, hairline + 2 * pxWing, hairline);
        g.fillRect(x - pxWing, y + h - hairline, hairline + 2 * pxWing, hairline);

        // draw
        double q = getInsertionQuality(alignment, block, renderOptions);
        if ( !Double.isNaN(q) && renderOptions.isIndelQualColoring() ) {
            Color currentColor = g.getColor();
            g.setColor(new Color(indelColorMap.getColor((int) q)));
            g.fillRect(x - pxWing, (int) (y + (h - hairline) * ((42 - q) / 42)) - 1, hairline + 2 * pxWing, hairline * 2);
            g.setColor(currentColor);
        }
    }

    public void renderSmallInsertionWings(Alignment alignment,
                                          AlignmentBlock block,
                                          RenderContext context,
                                          int pxH, int pxTop, int pxRight, int pxLeft,
                                          AlignmentTrack.RenderOptions renderOptions) {

        int pxWing = (pxH > 10 ? 2 : 1);  // width of the cursor "wing"
        Graphics2D g = context.getGraphics();

        // adjust wing and hairline
        int hairline = 2;
        double locScale = context.getScale();
        if ( renderOptions.isIndelQualColoring() ) {
            pxWing = Math.min(pxWing, Math.max(1, (int) (1 / locScale)));
            hairline = Math.min(hairline, pxWing);
        }

        Color currentColor = g.getColor();
        g.setColor(AlignmentRenderer.purple);
        g.fillRect(pxLeft - pxWing, pxTop, pxRight - pxLeft + hairline * pxWing, hairline);
        g.fillRect(pxLeft - pxWing, pxTop + pxH - hairline, pxRight - pxLeft + hairline * pxWing, hairline);
        g.setColor(currentColor);

        // draw
        double q = getInsertionQuality(alignment, block, renderOptions);
        if ( !Double.isNaN(q) && renderOptions.isIndelQualColoring() ) {
            currentColor = g.getColor();
            g.setColor(new Color(indelColorMap.getColor((int) q)));
            g.fillRect(pxLeft - pxWing, (int) (pxTop + (pxH - hairline) * ((42 - q) / 42)), pxRight - pxLeft + hairline * pxWing, hairline);
            g.setColor(currentColor);
        }
    }

    public void renderDeletionGap(Alignment alignment,
                                  Gap gap,
                                  int y, int h, int x, int w,
                                  RenderContext context,
                                  AlignmentTrack.RenderOptions renderOptions) {

        if ( !renderOptions.isIndelQualColoring() )
            return;

        // establish alignment blocks wrapping this gap (before and after)
        AlignmentBlock[] blocks = getGapWrappingBlocks(alignment, gap);
        if ( blocks == null )
            return;
        AlignmentBlock abPrev = blocks[0];
        AlignmentBlock abNext = blocks[1];

        // establish the nature of the gap
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        char        nextBlcokFirstBase = Character.toUpperCase((char)abNext.getBases().getByte(0));
        char        gapLastBase = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases() - 1));
        boolean     isForwardHmer = gapLastBase == nextBlcokFirstBase;
        char        prevBlcokLastBase = Character.toUpperCase((char)abPrev.getBases().getByte(abPrev.getBasesLength() - 1));
        char        gapFirstBase = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart()));
        boolean     isBackwardsHmer = gapFirstBase == prevBlcokLastBase;

        // we call something an hmer in this context only if it is hmer on exactly one side
        boolean     isHmer = isBackwardsHmer ^ isForwardHmer;

        // in addition, the hmer must span the whole gap
        if ( isHmer ) {
            isHmer = gapIsAllSameBase(genome, alignment, gap);
        }

        // figure out quality to plot - if any
        double q = Double.NaN;
        if ( isHmer ) {
            final SAMRecord record = ((SAMAlignment)alignment).getRecord();
            Hmer hmer = isForwardHmer
                    ? findHmer(record, abNext.getIndexOnRead(), (byte) gapLastBase, false, true)
                    : findHmer(record, abPrev.getIndexOnRead() + abPrev.getBasesLength() - 1, (byte) gapFirstBase, true, false);

            if ( hmer.size() >= getMC(record) ) {
                // HMER - length is at least max-hmer
                q = MIN_POSSIBLE_QUALITY;
            } else {
                // HMER - otherwise try TP
                q = getQualityFromTP(record, hmer, gap.getnBases(), renderOptions);
            }
        } else {
            if ( gap.getnBases() == 1 ) {
                // NON-HMER, single base, try T0
                final SAMRecord record = ((SAMAlignment)alignment).getRecord();
                double qBefore = getQualityFromT0(record, abPrev, false);
                double qAfter = getQualityFromT0(record, abNext, true);
                if ( Double.isNaN(qBefore) )
                    q = qAfter;
                else if ( Double.isNaN(qAfter) )
                    q = qBefore;
                else
                    q = Math.max(qBefore, qAfter);
            }
        }
        if ( Double.isNaN(q) )
            return;

        // draw the marker
        int hairline = Math.min(2, (int) (1 / context.getScale()));
        Graphics2D g = context.getGraphics();
        Color currentColor = g.getColor();
        g.setColor(new Color(indelColorMap.getColor((int)q)));
        g.fillRect(x + (int) (w * q / 42) - hairline, y, hairline * 2, h);
        g.setColor(currentColor);
    }

    private double getInsertionQuality(Alignment alignment, AlignmentBlock block, AlignmentTrack.RenderOptions renderOptions) {

        // establish alignment blocks wrapping this insertion block (before and after)
        AlignmentBlock[] blocks = getBlockWrappingBlocks(alignment, block);
        if ( blocks == null )
            return Double.NaN;
        AlignmentBlock abPrev = blocks[0];
        AlignmentBlock abNext = blocks[1];

        // establish the nature of the block
        char        nextBlcokFirstBase = Character.toUpperCase((char)abNext.getBases().getByte(0));
        char        blockLastBase = Character.toUpperCase((char)block.getBases().getByte(block.getBasesLength() - 1));
        boolean     isForwardHmer = blockLastBase == nextBlcokFirstBase;
        char        prevBlcokLastBase = Character.toUpperCase((char)abPrev.getBases().getByte(abPrev.getBasesLength() - 1));
        char        blockFirstBase = Character.toUpperCase((char)block.getBases().getByte(0));
        boolean     isBackwardsHmer = blockFirstBase == prevBlcokLastBase;

        // we call something an hmer in this context only if it is hmer on exactly one side
        boolean     isHmer = isBackwardsHmer ^ isForwardHmer;

        // in addition, the hmer must span the whole gap
        if ( isHmer ) {
            isHmer = blockIsAllSameBase(block);
        }

        // if both side are non-hmer and size is 1 - special case
        SAMRecord record = ((SAMAlignment)alignment).getRecord();
        Double q = Double.NaN;
        if ( !isBackwardsHmer && !isForwardHmer && block.getBasesLength() == 1 ) {

            // treat insertion as an hmer on to itself
            Hmer hmer = new Hmer();
            hmer.start = hmer.end = block.getIndexOnRead();
            q = getQualityFromTP(record, hmer, -hmer.size(), renderOptions);
            return q;
        }

        // if not an hmer, nothing to print
        if ( !isHmer )
            return Double.NaN;

        if ( isForwardHmer ) {
            Hmer hmer = findHmer(record, abNext.getIndexOnRead(), (byte)nextBlcokFirstBase, true, true);
            q = getQualityFromTP(record, hmer, -hmer.backwardsSize, renderOptions);
        } else if ( isBackwardsHmer ) {
            Hmer hmer = findHmer(record, abPrev.getIndexOnRead() + abPrev.getBasesLength() - 1 , (byte)prevBlcokLastBase, true, true);
            q = getQualityFromTP(record, hmer, -hmer.forwardSize, renderOptions);
        }

        return q;
    }

    private AlignmentBlock[] getBlockWrappingBlocks(Alignment alignment, AlignmentBlock block) {

        AlignmentBlock abPrev = null;
        AlignmentBlock abNext = null;
        for ( final AlignmentBlock b : alignment.getAlignmentBlocks() ) {
            if ( b.getEnd()  == block.getStart() )
                abPrev = b;
            else if ( b.getStart() == block.getStart() )
                abNext = b;
        }

        if ( abPrev != null && abNext != null ) {
            return new AlignmentBlock[] {abPrev, abNext};
        } else {
            return null;
        }
    }

    private AlignmentBlock[] getGapWrappingBlocks(Alignment alignment, Gap gap) {

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
            return new AlignmentBlock[] {alignment.getAlignmentBlocks()[blockIndex], alignment.getAlignmentBlocks()[blockIndex + 1]};
        } else {
            return null;
        }
    }

    private double getQualityFromTP(SAMRecord record, Hmer hmer, int tpValue, AlignmentTrack.RenderOptions renderOptions) {

        // get quals and tp
        byte[]      tp = record.getByteArrayAttribute(ATTR_TP);

        // scan for tpValue
        if ( !renderOptions.isIndelQualUsesMin() ) {
            for (int ofs = hmer.start; ofs <= hmer.end; ofs++) {
                if (tp[ofs] == tpValue) {
                    // all non-central entries must be doubled due to the symmetric nature of the quality string
                    final double q = record.getBaseQualities()[ofs];
                    final boolean isCenteralBaseInHmer = (ofs - hmer.start) == (hmer.end - ofs);
                    if (isCenteralBaseInHmer) {
                        return q;
                    } else {
                        return -10.0 * Math.log10(Math.pow(10, q / -10.0) * 2);
                    }
                }
            }

            // if here, not found. return highest quality value found on read
            final byte[] quals = record.getBaseQualities();
            if (quals.length == 0)
                return Double.NaN;
            double q = Double.MIN_VALUE;
            for (byte qual : quals) {
                q = Math.max(qual, q);
            }
            return q;
        } else {
            double q = Double.MAX_VALUE;
            for (int ofs = hmer.start; ofs <= hmer.end; ofs++) {
                q = Math.min(q, record.getBaseQualities()[ofs]);
            }
            return q;
        }
    }

    private double getQualityFromT0(SAMRecord record, AlignmentBlock block, boolean delIsBeforeBlock) {
        final String t0Attribute = record.getStringAttribute(TAG_T0);
        if ( t0Attribute == null ) {
            return Double.NaN;
        } else {
            final byte[] t0 = t0Attribute.getBytes();
            final int t0Index = delIsBeforeBlock ? block.getIndexOnRead() : (block.getIndexOnRead() + block.getLength() - 1);
            if ( t0Index < 0 || t0Index >= t0.length ) {
                return Double.NaN;
            } else {
                return t0[t0Index] - 33;
            }
        }
    }

    private int getMC(SAMRecord record) {
        try {
            if ( record.getReadGroup() == null )
                return 0;
            return Integer.parseInt(record.getReadGroup().getAttribute("mc"));
        } catch (Exception e) {
            return 0;
        }
    }

    private Hmer findHmer(SAMRecord record, int start, byte base, boolean walkBackwards, boolean walkForward) {

        Hmer hmer = new Hmer();
        hmer.end = hmer.start = start;
        byte[] bases = record.getReadBases();
        if ( walkBackwards ) {
            while (hmer.start > 0 && bases[hmer.start - 1] == base) {
                hmer.start--;
                hmer.backwardsSize++;
            }
        }
        if ( walkForward ) {
            while ((hmer.end + 1) < bases.length && bases[hmer.end + 1] == base) {
                hmer.end++;
                hmer.forwardSize++;
            }
        }

        return hmer;
    }

    private boolean gapIsAllSameBase(Genome genome, Alignment alignment, Gap gap) {
        final byte base = genome.getReference(alignment.getChr(), gap.getStart());
        for ( int i = 1 ; i < gap.getnBases() ; i++ ) {
            if ( genome.getReference(alignment.getChr(), gap.getStart() + i) != base ) {
                return false;
            }
        }
        return true;
    }

    private boolean blockIsAllSameBase(AlignmentBlock block) {
        final byte base = block.getBases().getByte(0);
        for ( int i = 1 ; i < block.getBasesLength() ; i++ ) {
            if ( block.getBases().getByte(i) != base ) {
                return false;
            }
        }
        return true;
    }


}