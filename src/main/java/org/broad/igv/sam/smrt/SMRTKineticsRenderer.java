package org.broad.igv.sam.smrt;

import org.broad.igv.sam.*;

import java.awt.*;

public class SMRTKineticsRenderer {
    public static void drawSmrtKinetics(Alignment alignment, double bpStart, double locScale, Rectangle rowRect, Graphics g, AlignmentTrack.ColorOption colorOption) {
        int dX;
        short[] smrtFrameCounts;
        SMRTKinetics smrtKinetics = alignment.getSmrtKinetics();
        if (AlignmentTrack.ColorOption.SMRT_SUBREAD_IPD == colorOption) {
            smrtFrameCounts = smrtKinetics.getSmrtSubreadIpd();
        } else if (AlignmentTrack.ColorOption.SMRT_SUBREAD_PW == colorOption) {
            smrtFrameCounts = smrtKinetics.getSmrtSubreadPw();
        } else if (AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD == colorOption || AlignmentTrack.ColorOption.SMRT_CCS_REV_IPD == colorOption) {
            final boolean isForwardStrand = (AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD == colorOption);
            smrtFrameCounts = smrtKinetics.getSmrtCcsIpd(isForwardStrand);
        } else {
            final boolean isForwardStrand = (AlignmentTrack.ColorOption.SMRT_CCS_FWD_PW == colorOption);
            smrtFrameCounts = smrtKinetics.getSmrtCcsPw(isForwardStrand);
        }

        if (smrtFrameCounts != null) {
            // Compute bounds
            int pY = (int) rowRect.getY();
            int dY = (int) rowRect.getHeight();
            dX = (int) Math.max(1, (1.0 / locScale));

            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {
                ByteSubarray bases = block.getBases();
                final int startOffset = bases.startOffset;
                final int stopOffset = startOffset + bases.length;
                for (int i = startOffset; i < stopOffset; i++) {
                    g.setColor(getSmrtFrameCountColor(smrtFrameCounts[i]));

                    int blockIdx = i - block.getBases().startOffset;
                    int pX = (int) ((block.getStart() + blockIdx - bpStart) / locScale);

                    // Don't draw out of clipping rect
                    if (pX > rowRect.getMaxX()) {
                        break;
                    } else if (pX + dX < rowRect.getX()) {
                        continue;
                    }

                    // Expand narrow width to make more visible
                    if (dX < 3) {
                        dX = 3;
                        pX--;
                    }
                    g.fillRect(pX, pY, dX, Math.max(1, dY - 2));
                }
            }
        }
    }

    /**
     * Set base color for SMRT sequencing frame count
     * <p>
     * This color scheme uses multiple transitions to help visually distinguish a large range of time intervals.
     * <p>
     * Color scheme:
     * The color starts out as transparent at a frame count of 0.
     * - Color transition 1 goes from transparent to opaque red
     * - Color transition 2 goes from red to yellow
     * - Color transition 3 goes from yellow to cyan
     * <p>
     * Transition 3 will mostly be visible when uncompressed frame counts are used. It helps to visually distinguish
     * exceptional polymerase stalling events.
     */
    private static Color getSmrtFrameCountColor(short shortFrameCount) {
        final int transition1MaxFrameCount = 100;
        final int transition2MaxFrameCount = 600;
        final int transition3MaxFrameCount = 6000;
        final Color color1 = Color.red;
        final Color color2 = Color.yellow;
        final Color color3 = Color.cyan;

        final int frameCount = Short.toUnsignedInt(shortFrameCount);
        final int alpha = Math.min(255, (frameCount * 255) / transition1MaxFrameCount);
        Color blendedColor;
        if (frameCount <= transition1MaxFrameCount) {
            blendedColor = color1;
        } else if (frameCount <= transition2MaxFrameCount) {
            float color2Fraction = (frameCount - transition1MaxFrameCount) / ((float) (transition2MaxFrameCount - transition1MaxFrameCount));
            blendedColor = blendColors(color1, color2, color2Fraction);
        } else if (frameCount <= transition3MaxFrameCount) {
            float color3Fraction = (frameCount - transition2MaxFrameCount) / ((float) (transition3MaxFrameCount - transition2MaxFrameCount));
            blendedColor = blendColors(color2, color3, color3Fraction);
        } else {
            blendedColor = color3;
        }
        return new Color(blendedColor.getRed(), blendedColor.getGreen(), blendedColor.getBlue(), alpha);
    }


    private static int blendColorValue(int v1, int v2, float v2Frac) {
        return Math.min(255, (int) (v1 + (v2 - v1) * v2Frac));
    }

    private static Color blendColors(Color c1, Color c2, float c2Frac) {
        assert (c2Frac >= 0.0 && c2Frac <= 1.0);
        return new Color(blendColorValue(c1.getRed(), c2.getRed(), c2Frac),
                blendColorValue(c1.getGreen(), c2.getGreen(), c2Frac),
                blendColorValue(c1.getBlue(), c2.getBlue(), c2Frac));
    }

}
