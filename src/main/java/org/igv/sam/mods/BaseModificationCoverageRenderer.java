package org.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;
import org.igv.sam.AlignmentCounts;
import org.igv.sam.AlignmentTrack.ColorOption;
import org.igv.track.RenderContext;

import java.awt.*;
import java.util.*;
import java.util.List;

public class BaseModificationCoverageRenderer {

    public static void drawModifications(RenderContext context,
                                         int pX,
                                         int pBottom,
                                         int dX,
                                         int barHeight,
                                         int pos,
                                         AlignmentCounts alignmentCounts,
                                         ColorOption colorOption,
                                         BaseModficationFilter filter,
                                         float threshold,
                                         Set<String> simplexModifications) {


        BaseModificationCounts modificationCounts = alignmentCounts.getModifiedBaseCounts();

        if (modificationCounts != null) {

            Set<BaseModificationKey> allModificationKeys = modificationCounts.getAllModificationKeys();
            List<BaseModificationKey> sortedKeys = new ArrayList<>(allModificationKeys);
            Collections.sort(sortedKeys);

            // Color bar modification counts (# of modifications > threshold)
            int total = alignmentCounts.getTotalCount(pos);

            for (BaseModificationKey key : sortedKeys) {

                try {
                    if (filter != null && !filter.pass(key.modification, key.getCanonicalBase())) continue;
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
                if (key.modification.startsWith("NONE_") && colorOption != ColorOption.BASE_MODIFICATION_2COLOR)
                    continue;

                byte base = (byte) key.getBase();
                final byte compl = SequenceUtil.complement(base);
                int modifiable = base == 'N' ? total : alignmentCounts.getCount(pos, base) + alignmentCounts.getCount(pos, compl);
                int detectable = base == 'N' ? total : simplexModifications.contains(key.modification) ?
                        alignmentCounts.getPosCount(pos, base) + alignmentCounts.getNegCount(pos, compl) :
                        modifiable;

                if (detectable == 0) continue;  //No informative reads

                int count = modificationCounts.getCount(pos, key, threshold, colorOption == ColorOption.BASE_MODIFICATION_2COLOR);
                if (count == 0) continue;

                float modFraction = (((float) modifiable) / total) * (((float) count) / detectable);
                int modHeight = Math.round(modFraction * barHeight);

                int likelihoodSum = modificationCounts.getLikelihoodSum(pos, key, threshold, colorOption == ColorOption.BASE_MODIFICATION_2COLOR);
                int averageLikelihood = (int) ((double) likelihoodSum) / count;

                int baseY = pBottom - modHeight;
                Color modColor = BaseModificationColors.getModColor(key.modification, averageLikelihood, colorOption);
                Graphics2D graphics = context.getGraphics();
                graphics.setColor(modColor);
                graphics.fillRect(pX, baseY, dX, modHeight);
                pBottom = baseY;

            }
        }
    }

    static Map<String, Integer> modRankOrder;


}
