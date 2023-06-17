package org.broad.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.AlignmentTrack.ColorOption;
import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class BaseModificationCoverageRenderer {

    public static void drawModifications(RenderContext context,
                                         int pX,
                                         int pBottom,
                                         int dX,
                                         int barHeight,
                                         int pos,
                                         AlignmentCounts alignmentCounts,
                                         ColorOption colorOption) {

        switch (colorOption) {
            case BASE_MODIFICATION_5MC:
                draw5MC(context, pX, pBottom, dX, barHeight, pos, alignmentCounts, false);
                break;
            case BASE_MODIFICATION_C:
                draw5MC(context, pX, pBottom, dX, barHeight, pos, alignmentCounts, true);
                break;
            case BASE_MODIFICATION_6MA:
                draw(context, pX, pBottom, dX, barHeight, pos, alignmentCounts, true);
                break;
            default:
                draw(context, pX, pBottom, dX, barHeight, pos, alignmentCounts, false);
        }
    }


    private static void draw(RenderContext context,
                             int pX,
                             int pBottom,
                             int dX,
                             int barHeight,
                             int pos,
                             AlignmentCounts alignmentCounts,
                             boolean onlyDraw6mA) {

        BaseModificationCounts modificationCounts = alignmentCounts.getModifiedBaseCounts();

        if (modificationCounts != null) {

            Graphics2D graphics = context.getGraphics();

            for (BaseModificationCounts.Key key : modificationCounts.getAllModifications()) {

                String modification = key.getModification();
                if (onlyDraw6mA) {
                    if (key.getCanonicalBase() != 'A' && key.getCanonicalBase() != 'T') continue;
                    if (!modification.equals("a")) continue;
                }

                // The number of modification calls, some of which might have likelihood of zero
                int modificationCount = modificationCounts.getCount(pos, key);

                if (barHeight > 0 && modificationCount > 0) {

                    byte base = (byte) key.getBase();
                    byte complement = SequenceUtil.complement(base);

                    // Count of bases at this location that could potentially be modified, accounting for strand
                    int baseCount = alignmentCounts.getPosCount(pos, base) + alignmentCounts.getNegCount(pos, complement);

                    int calledBarHeight = (int) ((((float) modificationCount) / baseCount) * barHeight);
                    Color modColor = BaseModificationColors.getModColor(modification, (byte) 255, ColorOption.BASE_MODIFICATION);

                    float averageLikelihood = (float) (modificationCounts.getLikelhoodSum(pos, key)) / (modificationCount * 255);
                    int modHeight = (int) (averageLikelihood * calledBarHeight);

                    // Generic modification
                    float threshold = PreferencesManager.getPreferences().getAsFloat("SAM.BASEMOD_THRESHOLD");
                    if (averageLikelihood > threshold && modHeight > 0) {
                        int baseY = pBottom - modHeight;
                        graphics.setColor(modColor);
                        graphics.fillRect(pX, baseY, dX, modHeight);
                        pBottom = baseY;
                    }
                }
            }
        }
    }

    private static void draw5MC(RenderContext context,
                                int pX,
                                int pBottom,
                                int dX,
                                int barHeight,
                                int pos,
                                AlignmentCounts alignmentCounts,
                                boolean allMods) {

        BaseModificationCounts modificationCounts = alignmentCounts.getModifiedBaseCounts();

        if (modificationCounts != null) {

            final byte base = (byte) 'C';
            final byte complement = (byte) 'G';

            Map<String, Integer> likelihoodSums = new HashMap<>();
            Map<String, Integer> modCounts = new HashMap<>();
            for (BaseModificationCounts.Key key : modificationCounts.getAllModifications()) {

                // This coloring mode is exclusively for "C" modifications
                if (key.getCanonicalBase() != 'C') continue;
                if (key.getModification().equals("m") || allMods) {
                    String mod = key.getModification();
                    final int count = modificationCounts.getCount(pos, key);
                    if (count > 0) {
                        modCounts.put(mod, count);
                        final int likelhoodSum = modificationCounts.getLikelhoodSum(pos, key);
                        likelihoodSums.put(mod, likelhoodSum);
                    }
                }
            }

            if (likelihoodSums.size() > 0) {

                // Count of bases at this location that could potentially be modified
                double modifiableBaseCount = alignmentCounts.getPosCount(pos, base) + alignmentCounts.getNegCount(pos, complement);

                // Compute "snp factor", ratio of count of base calls that could be modfied (on either strand) to
                // total count. This is normally close to 1, but can be less due non CG bases at this location (e.g. snps)
                double cgCount = alignmentCounts.getCount(pos, base) + alignmentCounts.getCount(pos, complement);
                double snpFactor = cgCount / alignmentCounts.getTotalCount(pos);

                double calledBarHeight = snpFactor * barHeight;
                double t = modifiableBaseCount * 255;   // If all bases are called this is the total sum of all likelihoods, including "no mod" likelihood

                // Likelihood of no modification.  It is assumed that the likelihood of no modification == (1 - sum(likelihood))
                int c = Collections.max(modCounts.values());
                double noModProb = c * 255;
                for (String m : modCounts.keySet()) {
                    if (likelihoodSums.containsKey(m)) {
                        noModProb -= likelihoodSums.get(m);
                    }
                }

                // Draw "no mod" bar
                int noModHeight = (int) Math.round((noModProb / t) * calledBarHeight);
                int baseY = pBottom - noModHeight;
                Graphics2D graphics = context.getGraphics();
                graphics.setColor(BaseModificationColors.noModColor5MC);
                graphics.fillRect(pX, baseY, dX, noModHeight);

                // Loop through modifications drawing bar for each
                String[] orderedMods = likelihoodSums.keySet().toArray(new String[0]);
                Arrays.sort(orderedMods, (o1, o2) -> -1 * o1.compareTo(o2));
                for (String m : orderedMods) {
                    Color mColor = BaseModificationColors.getModColor(m, (byte) 255, ColorOption.BASE_MODIFICATION_5MC);
                    int mModHeight = (int) Math.round(((likelihoodSums.get(m)) / t) * calledBarHeight);
                    if (mModHeight > 0) {
                        baseY -= mModHeight;
                        graphics.setColor(mColor);
                        graphics.fillRect(pX, baseY, dX, mModHeight);
                    }
                }
            }
        }
    }
}
