package org.broad.igv.sam.mods;

import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.AlignmentTrack.ColorOption;
import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.*;

public class BaseModificationCoverageRenderer {



    public static void drawModifications(RenderContext context,
                                         int pX,
                                         int pBottom,
                                         int dX,
                                         int barHeight,
                                         int pos,
                                         AlignmentCounts alignmentCounts,
                                         ColorOption colorOption,
                                         String basemodFilter) {

        switch (colorOption) {
            case BASE_MODIFICATION_C:
                draw5MC(context, pX, pBottom, dX, barHeight, pos, alignmentCounts, basemodFilter);
                break;
            default:
                draw(context, pX, pBottom, dX, barHeight, pos, alignmentCounts, colorOption, basemodFilter);
        }
    }

    /*
    chr11:119,094,758
Total count: 26
A : 0
C : 26 (100%, 16+, 10- )
G : 0
T : 0
N : 0Modification: m (C+, 15 @ 47%)

     */

    private static void draw(RenderContext context,
                             int pX,
                             int pBottom,
                             int dX,
                             int barHeight,
                             int pos,
                             AlignmentCounts alignmentCounts,
                             ColorOption colorOption,
                             String filter) {

        boolean debug = pos == 119094762;

        BaseModificationCounts modificationCounts = alignmentCounts.getModifiedBaseCounts();

        if (modificationCounts != null) {

            final  BaseModificationKey mCKey =  BaseModificationKey.getKey('C', '+', "m");
            final  BaseModificationKey mGKey =  BaseModificationKey.getKey('G', '-', "m");

            Graphics2D graphics = context.getGraphics();

            Set<BaseModificationKey> allModificationKeys = modificationCounts.getAllModificationKeys();
            boolean cpgMode = allModificationKeys.contains(mCKey) && !allModificationKeys.contains(mGKey);

            // Merge complementary sets (same modification, opposite strands)
            Map<String, Float> likelihoodSums = new LinkedHashMap<>();
            for (BaseModificationKey key : allModificationKeys) {
                String modification = key.getModification();
                if(filter != null && !filter.equals(modification)) continue;
                float currentCount = likelihoodSums.containsKey(modification) ? likelihoodSums.get(modification) : 0;
                likelihoodSums.put(modification, currentCount + modificationCounts.getLikelhoodSum(pos, key));
            }

            // Color bar by likelihood weighted count
            int total = alignmentCounts.getTotalCount(pos);

            int sumModHeight = 0;
            for (String modification : likelihoodSums.keySet()) {

                if (likelihoodSums.get(modification) > 0) {

                    float modFraction;
                    if (cpgMode & "m".equals(modification)) {
                        // Special mode for out-of-spec 5mC CpG convention.  Calls are recorded for C+ only.
                        // We adjust the height of the bar to account for the missing G- calls.
                        final byte base = (byte) 'C';
                        final byte compl = (byte) 'G';
                        final int posCount = alignmentCounts.getPosCount(pos, base);
                        final int negCount = alignmentCounts.getNegCount(pos, compl);

                        // Is this a "C" or "G" reference site?  We want to determine this from the counts data
                        // directly, this should work for all but pathological edge cases
                        final int cCounts = alignmentCounts.getCount(pos, base);
                        final int gCounts = alignmentCounts.getCount(pos, compl);
                        final boolean cSite = cCounts > gCounts;
                        final int referenceCounts = cSite ? cCounts : gCounts;
                        final int strandCount = cSite ? posCount : negCount;
                        modFraction = (((float) referenceCounts) / total) * (likelihoodSums.get(modification) / (strandCount * 255f));
                    } else {
                        modFraction = likelihoodSums.get(modification) / (total * 255f);
                    }
                    int modHeight = Math.round(modFraction * barHeight);

                    int baseY = pBottom - modHeight;
                    Color modColor = BaseModificationColors.getModColor(modification, (byte) 255, colorOption);
                    graphics.setColor(modColor);
                    graphics.fillRect(pX, baseY, dX, modHeight);
                    pBottom = baseY;

                    if(debug) System.out.println("draw  " + modification + "  " + modHeight);
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
                                String basemodFilter) {

        boolean debug = pos == 119094768;

        BaseModificationCounts modificationCounts = alignmentCounts.getModifiedBaseCounts();

        if (modificationCounts != null) {

            final byte base = (byte) 'C';
            final byte compl = (byte) 'G';

            Map<String, Integer> likelihoodSums = new HashMap<>();
            Map<String, Integer> modCounts = new HashMap<>();
            for (BaseModificationKey key : modificationCounts.getAllModificationKeys()) {

                // This coloring mode is exclusively for "C" modifications
                if (key.getCanonicalBase() != 'C') continue;
                if (key.getModification().equals(basemodFilter) || basemodFilter == null) {
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

                // Special mode for out-of-spec 5mC CpG convention.
                // Calls are made for the CG dinucleotide and only recorded on 1 strand.  We adjust the height
                // of the bar to account for the missing G- calls.  This is an approximation and assumes the
                // distribution of calls is ~ equal on both strands.

                // Is this a "C" or "G" reference site?  We want to determine this from the counts data directly,
                // this should work for all but pathological edge cases
                final int cCounts = alignmentCounts.getCount(pos, base);
                final int gCounts = alignmentCounts.getCount(pos, compl);
                final boolean cSite = cCounts > gCounts;
                final int referenceCounts = cSite ? cCounts : gCounts;

               // Compute "snp factor", ratio of count of base calls that could be modified  to
                // total count. This is normally close to 1, but can be less due non CG bases at this location (e.g. snps)
                double snpFactor = ((double) referenceCounts) / alignmentCounts.getTotalCount(pos);

                double calledBarHeight = snpFactor * barHeight;
                double t = referenceCounts * 255;   // If all bases are called this is the total sum of all likelihoods, including "no mod" likelihood

                // Likelihood of no modification.  It is assumed that the likelihood of no modification == (1 - sum(likelihood))
                //int c = Collections.max(modCounts.values());
                double noModProb = t;
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
                    Color mColor = BaseModificationColors.getModColor(m, (byte) 255, ColorOption.BASE_MODIFICATION_C);
                    int mModHeight = (int) Math.round(((likelihoodSums.get(m)) / t) * calledBarHeight);
                    if (mModHeight > 0) {
                        baseY -= mModHeight;
                        graphics.setColor(mColor);
                        graphics.fillRect(pX, baseY, dX, mModHeight);
                    }
                    if(debug) System.out.println("draw5mC  " + m + "  " + mModHeight);
                }
            }
        }
    }
}
