package org.igv.sashimi;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class SashimiCoordinateMapTest {

    private static final double EXPONENT = SashimiCoordinateMap.DEFAULT_EXPONENT;

    @Test
    public void intersectJunctions() {
        // Exon skipping -- the long junction does not shrink the skipped exon at 400-600
        List<int[]> regions = SashimiCoordinateMap.intersectJunctions(List.of(
                new int[]{100, 1000}, new int[]{600, 1000}, new int[]{100, 400}));
        assertEquals(2, regions.size());
        assertArrayEquals(new int[]{100, 400}, regions.get(0));
        assertArrayEquals(new int[]{600, 1000}, regions.get(1));

        // Touching junctions do not overlap
        regions = SashimiCoordinateMap.intersectJunctions(List.of(new int[]{100, 200}, new int[]{200, 300}));
        assertEquals(2, regions.size());

        assertTrue(SashimiCoordinateMap.intersectJunctions(List.of()).isEmpty());
    }

    @Test
    public void toPlot() {
        SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(
                new int[]{1000, 11000}, new int[]{11100, 12100}), EXPONENT);
        double shrunk1 = Math.pow(10000, EXPONENT);
        double shrunk2 = Math.pow(1000, EXPONENT);

        assertEquals(500, map.toPlot(500), 1e-9);
        assertEquals(1000, map.toPlot(1000), 1e-9);
        assertEquals(1000 + shrunk1 / 2, map.toPlot(6000), 1e-9);
        assertEquals(1000 + shrunk1, map.toPlot(11000), 1e-9);
        // The exon between the introns keeps its width
        assertEquals(1000 + shrunk1 + 100, map.toPlot(11100), 1e-9);
        assertEquals(1000 + shrunk1 + 100 + shrunk2 + 50, map.toPlot(12150), 1e-9);
    }

    @Test
    public void exponent() {
        List<int[]> junctions = List.of(new int[]{1000, 11000});

        // A larger exponent compresses less
        assertEquals(1000 + Math.pow(10000, 0.3), SashimiCoordinateMap.fromJunctions(junctions, 0.3).toPlot(11000), 1e-9);
        assertEquals(1000 + Math.pow(10000, 0.9), SashimiCoordinateMap.fromJunctions(junctions, 0.9).toPlot(11000), 1e-9);

        // Same regions at different exponents must not compare equal, or the plot would skip the rebuild
        assertNotEquals(SashimiCoordinateMap.fromJunctions(junctions, 0.3),
                SashimiCoordinateMap.fromJunctions(junctions, 0.9));
    }

    @Test
    public void toGenomicInvertsToPlot() {
        for (double exponent : new double[]{0.3, EXPONENT, 0.9}) {
            SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(
                    new int[]{1000, 11000}, new int[]{11100, 12100}), exponent);
            for (double position : new double[]{0, 999.5, 1000, 4321, 10999, 11000, 11050, 11100, 12000, 20000}) {
                assertEquals(position, map.toGenomic(map.toPlot(position)), 1e-6);
            }
        }
    }

    @Test
    public void identity() {
        SashimiCoordinateMap map = SashimiCoordinateMap.identity();
        assertEquals(12345.5, map.toPlot(12345.5), 0);
        assertEquals(12345.5, map.toGenomic(12345.5), 0);
        assertEquals(map, SashimiCoordinateMap.fromJunctions(List.of(), EXPONENT));
    }

    @Test
    public void equality() {
        SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(new int[]{100, 400}), EXPONENT);
        assertEquals(map, SashimiCoordinateMap.fromJunctions(List.of(new int[]{100, 400}), EXPONENT));
        assertNotEquals(map, SashimiCoordinateMap.fromJunctions(List.of(new int[]{100, 500}), EXPONENT));
    }
}
