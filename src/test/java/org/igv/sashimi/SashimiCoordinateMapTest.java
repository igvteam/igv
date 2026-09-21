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
    public void drawnLength() {
        // An intron of length L is drawn with width L^exponent
        assertEquals(Math.pow(10000, EXPONENT), SashimiCoordinateMap.drawnLength(10000, EXPONENT), 1e-9);

        // An exponent of 1 draws introns at their true width
        for (double length : new double[]{99, 1500, 123000}) {
            assertEquals(length, SashimiCoordinateMap.drawnLength(length, 1), 1e-9);
        }

        // Longer introns are compressed harder, in proportion
        assertTrue(SashimiCoordinateMap.drawnLength(100000, EXPONENT) / 100000
                < SashimiCoordinateMap.drawnLength(1000, EXPONENT) / 1000);

        // Monotonic in both length and exponent
        assertTrue(SashimiCoordinateMap.drawnLength(50000, EXPONENT) > SashimiCoordinateMap.drawnLength(10000, EXPONENT));
        assertTrue(SashimiCoordinateMap.drawnLength(50000, 0.5) > SashimiCoordinateMap.drawnLength(50000, 0.3));
    }

    @Test
    public void toPlot() {
        // A long intron, a shorter one, and a 100 bp exon between them
        SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(
                new int[]{1000, 101000}, new int[]{101100, 102100}), EXPONENT);
        double shrunkLong = SashimiCoordinateMap.drawnLength(100000, EXPONENT);
        double shrunkShort = SashimiCoordinateMap.drawnLength(1000, EXPONENT);

        assertEquals(500, map.toPlot(500), 1e-9);
        assertEquals(1000, map.toPlot(1000), 1e-9);
        assertEquals(1000 + shrunkLong / 2, map.toPlot(51000), 1e-9);
        assertEquals(1000 + shrunkLong, map.toPlot(101000), 1e-9);
        // The exon between the introns keeps its width
        assertEquals(1000 + shrunkLong + 100, map.toPlot(101100), 1e-9);
        assertEquals(1000 + shrunkLong + 100 + shrunkShort, map.toPlot(102100), 1e-9);
    }

    @Test
    public void exponent() {
        List<int[]> junctions = List.of(new int[]{1000, 101000});

        // A larger exponent compresses less
        assertTrue(SashimiCoordinateMap.fromJunctions(junctions, 0.5).toPlot(101000)
                > SashimiCoordinateMap.fromJunctions(junctions, 0.3).toPlot(101000));

        // Same regions at different exponents must not compare equal, or the plot would skip the rebuild
        assertNotEquals(SashimiCoordinateMap.fromJunctions(junctions, 0.3),
                SashimiCoordinateMap.fromJunctions(junctions, 0.5));

        // An exponent of 1 leaves introns at their true width
        SashimiCoordinateMap uncompressed = SashimiCoordinateMap.fromJunctions(junctions, 1.0);
        for (double position : new double[]{0, 1000, 51000, 101000, 200000}) {
            assertEquals(position, uncompressed.toPlot(position), 1e-9);
        }
    }

    @Test
    public void toGenomicInvertsToPlot() {
        for (double exponent : new double[]{0, 0.3, 0.5, 1.0}) {
            SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(
                    new int[]{1000, 101000}, new int[]{101100, 102100}), exponent);
            for (double position : new double[]{0, 999.5, 1000, 4321, 100999, 101000, 101050, 101100, 102000, 200000}) {
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
}
