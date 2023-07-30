package org.broad.igv.sam.mods;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.AlignmentTrack;

import java.awt.*;
import java.util.List;

public class BaseModificationRenderer {

    public static void drawModifications(
            Alignment alignment,
            double bpStart,
            double locScale,
            Rectangle rowRect,
            Graphics g,
            AlignmentTrack.ColorOption colorOption,
            BaseModficationFilter filter) {

        List<BaseModificationSet> baseModificationSets = alignment.getBaseModificationSets();

        if (baseModificationSets != null) {
            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {
                // Compute bounds
                int pY = (int) rowRect.getY();
                int dY = (int) rowRect.getHeight();
                int dX = (int) Math.max(1, (1.0 / locScale));

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
                    boolean modificationFound = false;
                    for (BaseModificationSet bmSet : baseModificationSets) {
                        if (bmSet.containsPosition(i)) {
                            int lh = Byte.toUnsignedInt(bmSet.getLikelihoods().get(i));
                            // TODO ELSE if flag == '.' lh = 0
                            noModLh -= lh;
                            if (modification == null || lh > maxLh) {
                                modification = bmSet.getModification();
                                canonicalBase = bmSet.getCanonicalBase();
                                maxLh = lh;
                                modificationFound = modificationFound || filter == null || filter.pass(modification, canonicalBase);
                            }
                        }
                    }

                    if(modificationFound) {
                        Color c = null;
                        if (noModLh > 127 && colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR) {
                            c = BaseModificationColors.getModColor("NONE_" + canonicalBase, noModLh, colorOption);
                        } else if (modification != null && maxLh > 127 &&
                                (filter == null || filter.pass(modification, canonicalBase)))         /* || flag == . */ {
                            c = BaseModificationColors.getModColor(modification, maxLh, colorOption);
                        }
                        if (c != null) {
                            g.setColor(c);

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
        }
    }
}


