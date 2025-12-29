package org.igv.sam.mods;

import org.igv.sam.Alignment;
import org.igv.sam.AlignmentBlock;
import org.igv.sam.AlignmentTrack;
import org.igv.sam.InsertionMarker;

import java.awt.*;
import java.util.List;

public class BaseModificationRenderer {

    public static void drawModifications(
            Alignment alignment,
            double bpStart,
            double locScale,
            Rectangle rowRect,
            Graphics g,
            AlignmentTrack.RenderOptions renderOptions) {


        List<BaseModificationSet> baseModificationSets = alignment.getBaseModificationSets();

        if (baseModificationSets != null) {
            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {
                // Compute bounds
                drawBlock(bpStart, locScale, rowRect, g, renderOptions, baseModificationSets, block, alignment.isNegativeStrand());
            }
        }
    }

    public static void drawBlock(double bpStart, double locScale, Rectangle rowRect, Graphics g, AlignmentTrack.RenderOptions renderOptions, List<BaseModificationSet> baseModificationSets, AlignmentBlock block, boolean isNegativeStrand) {

        if(baseModificationSets == null) return;

        AlignmentTrack.ColorOption colorOption = renderOptions.getColorOption();
        BaseModficationFilter filter = renderOptions.getBasemodFilter();
        float threshold = renderOptions.getBasemodThreshold();
        boolean distinguishStrands = renderOptions.getBasemodDistinguishStrands();

        // Determine if we should leave a margin (same logic as AlignmentRenderer)
        boolean leaveMargin = renderOptions.getDisplayMode() != org.igv.track.Track.DisplayMode.SQUISHED;

        int pY = (int) rowRect.getY();
        int dY = (int) rowRect.getHeight();
        int dX = (int) Math.max(1, (1.0 / locScale));

        char alStrand = (isNegativeStrand ? '-' : '+');

        for (int i = block.getBases().startOffset; i < block.getBases().startOffset + block.getBases().length; i++) {

            int blockIdx = i - block.getBases().startOffset;
            int pX = (int) ((block.getStart() + blockIdx - bpStart) / locScale);

            // Don't draw out of clipping rect
            if (pX > rowRect.getMaxX()) {
                break;
            } else if (pX + dX < rowRect.getX()) {
                continue;
            }

            // Search all sets for modifications of this base.  Pick modification with largest likelhiood
            int maxLh = 0;
            int noModLh = 255;
            String modification = null;
            char canonicalBase = 0;
            char modStrand = 0;
            for (BaseModificationSet bmSet : baseModificationSets) {
                if (bmSet.containsPosition(i)) {
                    int lh = Byte.toUnsignedInt(bmSet.getLikelihoods().get(i));
                    noModLh -= lh;
                    if ((filter == null || filter.pass(bmSet.getModification(), canonicalBase)) && (modification == null || lh > maxLh)) {
                        modification = bmSet.getModification();
                        canonicalBase = bmSet.getCanonicalBase();
                        maxLh = lh;
                        modStrand = bmSet.getStrand();
                    }
                }
            }

            if (modification != null) {

                Color c = null;
                final float scaledThreshold = threshold * 255;
                if (noModLh > maxLh && colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR && noModLh >= scaledThreshold) {
                    c = BaseModificationColors.getModColor("NONE_" + canonicalBase, noModLh, colorOption);
                } else if (maxLh >= scaledThreshold) {
                    c = BaseModificationColors.getModColor(modification, maxLh, colorOption);
                }
                if (c != null) {
                    g.setColor(c);

                    // Expand narrow width to make more visible
                    if (dX < 3) {
                        dX = 3;
                        pX--;
                    }

                    int rectY = pY;
                    int rectHeight = Math.max(1, dY - (leaveMargin ? 2 : 0));

                    // Draw strand-specific half-height rectangles
                    // if "distinguish strands" setting is active
                    // Forward reference strand: upper half
                    // Reverse reference strand: lower half
                    // Only split by strand if height is sufficient (at least 4 pixels)
                    if (rectHeight >= 4 && modStrand != 0 && distinguishStrands) {
                        int halfHeight = rectHeight / 2;
                        if (modStrand == alStrand) {
                            // Forward reference strand: upper half
                            rectHeight = halfHeight;
                        } else {
                            // Reverse reference strand: lower half
                            rectY = pY + rectHeight - halfHeight;
                            rectHeight = halfHeight;
                        }
                    }

                    g.fillRect(pX, rectY, dX, rectHeight);
                }

            }
        }
    }
}


