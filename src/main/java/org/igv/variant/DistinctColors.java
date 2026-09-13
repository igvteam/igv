package org.igv.variant;

import org.igv.ui.color.ColorPalette;
import org.igv.ui.color.ColorUtilities;

import java.awt.Color;
import java.util.Collection;

/**
 * Picks a color for a new attribute value that can be told apart from the colors already in use -- by a track's
 * palette assignments, or by a scheme being edited.  Shared so both follow the same rules.
 */
final class DistinctColors {

    /**
     * Palette tried first.  Matches the default used by igv.js.
     */
    static final String PALETTE = "Set 1";

    /**
     * How far apart, in RGB, two colors have to be to count as distinguishable.
     */
    static final int MIN_DISTANCE = 60;

    /**
     * How many generated colors to try before settling for the furthest one found.  Also, therefore, how many
     * distinct colors can be in use before one is repeated.
     */
    static final int MAX_ATTEMPTS = 2000;

    private static final float[] SATURATION_LEVELS = {0.85f, 0.55f, 1.0f};
    private static final float[] BRIGHTNESS_LEVELS = {0.85f, 0.65f, 1.0f};

    private DistinctColors() {
    }

    /**
     * A color at least MIN_DISTANCE from everything in use.  Palette colors are preferred; once they are used up
     * or rejected, a non-repeating sequence is walked.  Policy on exhaustion: the first candidate at the required
     * distance is returned; failing that, the furthest; a color already in use is returned only if every one of
     * the MAX_ATTEMPTS candidates is -- which needs more distinct colors in use than there are candidates, since
     * the candidates are pairwise distinct.
     */
    static Color next(Collection<Color> used) {

        ColorPalette palette = ColorUtilities.getPalette(PALETTE);
        if (palette != null) {
            for (Color candidate : palette.getColors()) {
                if (minDistance(candidate, used) >= MIN_DISTANCE) {
                    return candidate;
                }
            }
        }

        Color best = null;
        double bestDistance = -1;
        for (int i = 0; i < MAX_ATTEMPTS; i++) {
            Color candidate = generated(i);
            double distance = minDistance(candidate, used);
            if (distance >= MIN_DISTANCE) {
                return candidate;
            }
            if (distance > bestDistance) {
                best = candidate;
                bestDistance = distance;
            }
        }
        return best;
    }

    /**
     * The i-th color of a sequence that does not repeat: the hue advances by the golden ratio each step, which
     * never returns to a previous hue, at a few saturation and brightness levels so consecutive candidates
     * differ in more than hue.  ({@link ColorUtilities#randomColor} is not usable here -- each channel is taken
     * modulo 215, so it has only 215 distinct colors, after which a search over it can only find duplicates.)
     */
    static Color generated(int i) {
        float hue = (float) ((i * 0.618033988749895) % 1.0);
        float saturation = SATURATION_LEVELS[i % SATURATION_LEVELS.length];
        float brightness = BRIGHTNESS_LEVELS[(i / SATURATION_LEVELS.length) % BRIGHTNESS_LEVELS.length];
        return Color.getHSBColor(hue, saturation, brightness);
    }

    /**
     * Distance from a color to the nearest of those in use, or Double.MAX_VALUE if none are.
     */
    static double minDistance(Color color, Collection<Color> used) {
        double min = Double.MAX_VALUE;
        for (Color c : used) {
            min = Math.min(min, distance(color, c));
        }
        return min;
    }

    static double distance(Color c1, Color c2) {
        int dr = c1.getRed() - c2.getRed();
        int dg = c1.getGreen() - c2.getGreen();
        int db = c1.getBlue() - c2.getBlue();
        return Math.sqrt(dr * dr + dg * dg + db * db);
    }
}
