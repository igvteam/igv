package org.igv.sashimi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Maps genomic positions to Sashimi plot positions, shrinking intronic regions.  Regions are chosen as in ggsashimi
 * (https://github.com/guigolab/ggsashimi):  overlapping junctions are intersected, and each intersection of length L
 * is drawn with length L^0.7.  Positions between regions keep their width.
 */
public class SashimiCoordinateMap {

    private static final double EXPONENT = 0.7;

    // Shrunk regions, sorted and non-overlapping, in genomic and plot coordinates
    private final double[] gStart;
    private final double[] gEnd;
    private final double[] pStart;
    private final double[] pEnd;

    public static SashimiCoordinateMap identity() {
        return new SashimiCoordinateMap(List.of());
    }

    /**
     * @param junctions junction {start, end} pairs
     */
    public static SashimiCoordinateMap fromJunctions(List<int[]> junctions) {
        return new SashimiCoordinateMap(intersectJunctions(junctions));
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

    SashimiCoordinateMap(List<int[]> regions) {
        int n = regions.size();
        gStart = new double[n];
        gEnd = new double[n];
        pStart = new double[n];
        pEnd = new double[n];
        double shift = 0;
        for (int i = 0; i < n; i++) {
            int[] r = regions.get(i);
            double length = r[1] - r[0];
            double shrunkLength = Math.pow(length, EXPONENT);
            gStart[i] = r[0];
            gEnd[i] = r[1];
            pStart[i] = r[0] - shift;
            pEnd[i] = pStart[i] + shrunkLength;
            shift += length - shrunkLength;
        }
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SashimiCoordinateMap)) return false;
        SashimiCoordinateMap other = (SashimiCoordinateMap) o;
        return Arrays.equals(gStart, other.gStart) && Arrays.equals(gEnd, other.gEnd);
    }

    @Override
    public int hashCode() {
        return 31 * Arrays.hashCode(gStart) + Arrays.hashCode(gEnd);
    }
}
