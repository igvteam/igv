package org.igv.sashimi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Maps genomic positions to Sashimi plot positions, shrinking intronic regions.  Regions are chosen as in ggsashimi
 * (https://github.com/guigolab/ggsashimi):  overlapping junctions are intersected.  Introns shorter than
 * {@link #THRESHOLD} are drawn at their true width, longer ones are compressed relative to the threshold, so that a
 * gene with one very long intron and many short ones keeps the short ones readable.  An exponent of 1 leaves all
 * introns at true width, an exponent of 0 draws every long intron at the threshold width.
 */
public class SashimiCoordinateMap {

    public static final double DEFAULT_EXPONENT = 0.3;

    /** Introns up to this length (bp) are never compressed */
    public static final double THRESHOLD = 2000;

    // Shrunk regions, sorted and non-overlapping, in genomic and plot coordinates
    private final double[] gStart;
    private final double[] gEnd;
    private final double[] pStart;
    private final double[] pEnd;

    public static SashimiCoordinateMap identity() {
        return new SashimiCoordinateMap(List.of(), DEFAULT_EXPONENT);
    }

    /**
     * @param junctions junction {start, end} pairs
     * @param exponent  a region of length L is drawn with width L^exponent
     */
    public static SashimiCoordinateMap fromJunctions(List<int[]> junctions, double exponent) {
        return new SashimiCoordinateMap(intersectJunctions(junctions), exponent);
    }

    /**
     * Port of ggsashimi's intersect_introns.  Junctions are sorted, and each run of overlapping junctions is reduced
     * to the region common to all of them.  Junctions that only touch do not overlap.
     */
    static List<int[]> intersectJunctions(List<int[]> junctions) {
        List<int[]> sorted = new ArrayList<>(junctions);
        sorted.sort(Comparator.<int[]>comparingInt(j -> j[0]).thenComparingInt(j -> j[1]));

        List<int[]> regions = new ArrayList<>();
        if (sorted.isEmpty()) {
            return regions;
        }
        int a = sorted.get(0)[0];
        int b = sorted.get(0)[1];
        for (int i = 1; i < sorted.size(); i++) {
            int c = sorted.get(i)[0];
            int d = sorted.get(i)[1];
            if (b > c) {
                b = Math.min(b, d);
                a = Math.max(a, c);
            } else {
                regions.add(new int[]{a, b});
                a = c;
                b = d;
            }
        }
        regions.add(new int[]{a, b});
        return regions;
    }

    SashimiCoordinateMap(List<int[]> regions, double exponent) {
        int n = regions.size();
        gStart = new double[n];
        gEnd = new double[n];
        pStart = new double[n];
        pEnd = new double[n];
        double shift = 0;
        for (int i = 0; i < n; i++) {
            int[] r = regions.get(i);
            double length = r[1] - r[0];
            double shrunkLength = drawnLength(length, exponent);
            gStart[i] = r[0];
            gEnd[i] = r[1];
            pStart[i] = r[0] - shift;
            pEnd[i] = pStart[i] + shrunkLength;
            shift += length - shrunkLength;
        }
    }

    /**
     * Width an intron of the given length is drawn with.  Introns at or below the threshold keep their width;
     * beyond it the excess is compressed, so drawn width grows as threshold * (length/threshold)^exponent.
     */
    static double drawnLength(double length, double exponent) {
        return length <= THRESHOLD ? length : THRESHOLD * Math.pow(length / THRESHOLD, exponent);
    }

    public double toPlot(double position) {
        int i = lastIndexAtOrBefore(gStart, position);
        if (i < 0) {
            return position;
        } else if (position >= gEnd[i]) {
            return pEnd[i] + (position - gEnd[i]);
        } else {
            return pStart[i] + (position - gStart[i]) * (pEnd[i] - pStart[i]) / (gEnd[i] - gStart[i]);
        }
    }

    public double toGenomic(double plotPosition) {
        int i = lastIndexAtOrBefore(pStart, plotPosition);
        if (i < 0) {
            return plotPosition;
        } else if (plotPosition >= pEnd[i]) {
            return gEnd[i] + (plotPosition - pEnd[i]);
        } else {
            return gStart[i] + (plotPosition - pStart[i]) * (gEnd[i] - gStart[i]) / (pEnd[i] - pStart[i]);
        }
    }

    private static int lastIndexAtOrBefore(double[] values, double x) {
        int i = Arrays.binarySearch(values, x);
        return i >= 0 ? i : -i - 2;
    }

    /**
     * Maps are equal if they place the same genomic regions at the same plot positions, so two maps computed from the
     * same junctions but with different exponents are not equal.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SashimiCoordinateMap)) return false;
        SashimiCoordinateMap other = (SashimiCoordinateMap) o;
        return Arrays.equals(gStart, other.gStart) && Arrays.equals(gEnd, other.gEnd)
                && Arrays.equals(pStart, other.pStart) && Arrays.equals(pEnd, other.pEnd);
    }

    @Override
    public int hashCode() {
        return 31 * (31 * Arrays.hashCode(gStart) + Arrays.hashCode(gEnd)) + Arrays.hashCode(pEnd);
    }
}
