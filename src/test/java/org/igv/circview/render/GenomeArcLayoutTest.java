package org.igv.circview.render;

import org.igv.circview.model.Assembly;
import org.igv.circview.model.Chromosome;
import org.junit.Test;

import java.awt.Color;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class GenomeArcLayoutTest {

    private static Assembly threeChrAssembly() {
        return new Assembly("test", "test", List.of(
                new Chromosome("1", 100, Color.RED),
                new Chromosome("2", 200, Color.GREEN),
                new Chromosome("3", 300, Color.BLUE)));
    }

    @Test
    public void spansPlusGapsCoverFullCircle() {
        double gapFraction = 0.02;
        GenomeArcLayout layout = new GenomeArcLayout(threeChrAssembly(), 0, 0, gapFraction);

        double totalSpan = 0;
        for (GenomeArcLayout.ChromosomeArc arc : layout.getArcs()) {
            totalSpan += arc.span();
        }
        double expectedAvailable = 2 * Math.PI * (1 - gapFraction);
        assertEquals(expectedAvailable, totalSpan, 1e-9);

        double totalGaps = layout.getGapAngle() * layout.getArcs().size();
        assertEquals(2 * Math.PI, totalSpan + totalGaps, 1e-9);
    }

    @Test
    public void spansAreProportionalToLength() {
        GenomeArcLayout layout = new GenomeArcLayout(threeChrAssembly(), 0, 0, 0.0);
        List<GenomeArcLayout.ChromosomeArc> arcs = layout.getArcs();
        // lengths 100:200:300 -> spans 1:2:3
        assertEquals(arcs.get(0).span() * 2, arcs.get(1).span(), 1e-9);
        assertEquals(arcs.get(0).span() * 3, arcs.get(2).span(), 1e-9);
    }

    @Test
    public void bpToAngleMatchesArcEndpoints() {
        GenomeArcLayout layout = new GenomeArcLayout(threeChrAssembly(), 0, 0, 0.02);
        GenomeArcLayout.ChromosomeArc arc2 = layout.getArcs().get(1); // chr "2", length 200

        assertEquals(arc2.startAngle, layout.bpToAngle("2", 0), 1e-9);
        assertEquals(arc2.endAngle, layout.bpToAngle("2", 200), 1e-9);
        double mid = (arc2.startAngle + arc2.endAngle) / 2;
        assertEquals(mid, layout.bpToAngle("2", 100), 1e-9);
    }

    @Test
    public void bpToAngleIsMonotonicWithinChromosome() {
        GenomeArcLayout layout = new GenomeArcLayout(threeChrAssembly(), 0, 0, 0.02);
        assertTrue(layout.bpToAngle("3", 10) < layout.bpToAngle("3", 250));
    }

    @Test
    public void unknownChromosome() {
        GenomeArcLayout layout = new GenomeArcLayout(threeChrAssembly(), 0, 0, 0.02);
        assertFalse(layout.hasChromosome("ZZ"));
        assertTrue(Double.isNaN(layout.bpToAngle("ZZ", 10)));
    }
}
