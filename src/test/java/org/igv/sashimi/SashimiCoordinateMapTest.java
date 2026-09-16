package org.igv.sashimi;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class SashimiCoordinateMapTest {

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
                new int[]{1000, 11000}, new int[]{11100, 12100}));
        double shrunk1 = Math.pow(10000, 0.7);
        double shrunk2 = Math.pow(1000, 0.7);

        assertEquals(500, map.toPlot(500), 1e-9);
        assertEquals(1000, map.toPlot(1000), 1e-9);
        assertEquals(1000 + shrunk1 / 2, map.toPlot(6000), 1e-9);
        assertEquals(1000 + shrunk1, map.toPlot(11000), 1e-9);
        // The exon between the introns keeps its width
        assertEquals(1000 + shrunk1 + 100, map.toPlot(11100), 1e-9);
        assertEquals(1000 + shrunk1 + 100 + shrunk2 + 50, map.toPlot(12150), 1e-9);
    }

    @Test
    public void toGenomicInvertsToPlot() {
        SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(
                new int[]{1000, 11000}, new int[]{11100, 12100}));
        for (double position : new double[]{0, 999.5, 1000, 4321, 10999, 11000, 11050, 11100, 12000, 20000}) {
            assertEquals(position, map.toGenomic(map.toPlot(position)), 1e-6);
        }
    }

    @Test
    public void identity() {
        SashimiCoordinateMap map = SashimiCoordinateMap.identity();
        assertEquals(12345.5, map.toPlot(12345.5), 0);
        assertEquals(12345.5, map.toGenomic(12345.5), 0);
        assertEquals(map, SashimiCoordinateMap.fromJunctions(List.of()));
    }

    @Test
    public void equality() {
        SashimiCoordinateMap map = SashimiCoordinateMap.fromJunctions(List.of(new int[]{100, 400}));
        assertEquals(map, SashimiCoordinateMap.fromJunctions(List.of(new int[]{100, 400})));
        assertNotEquals(map, SashimiCoordinateMap.fromJunctions(List.of(new int[]{100, 500})));
    }
}
