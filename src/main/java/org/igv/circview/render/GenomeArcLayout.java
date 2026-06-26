package org.igv.circview.render;

import org.igv.circview.model.Assembly;
import org.igv.circview.model.Chromosome;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Maps genomic coordinates to angles around the circle and angles to points.
 *
 * <p>Chromosomes are laid out clockwise starting at the top (12 o'clock), each
 * spanning an angle proportional to its length, with a small uniform gap between
 * neighbours. This class is pure geometry (no Swing) so it can be unit-tested.
 *
 * <p>Angle convention: an angle {@code a} maps to the point
 * {@code (cx + r*cos(a), cy + r*sin(a))}. Because screen y grows downward,
 * starting at {@code -PI/2} (top) and increasing the angle sweeps clockwise.
 */
public final class GenomeArcLayout {

    /** Angular extent of one chromosome's arc, in radians. */
    public static final class ChromosomeArc {
        public final Chromosome chromosome;
        public final double startAngle;
        public final double endAngle;

        ChromosomeArc(Chromosome chromosome, double startAngle, double endAngle) {
            this.chromosome = chromosome;
            this.startAngle = startAngle;
            this.endAngle = endAngle;
        }

        public double span() {
            return endAngle - startAngle;
        }
    }

    private static final double START_ANGLE = -Math.PI / 2.0; // top of circle

    private final double centerX;
    private final double centerY;
    private final List<ChromosomeArc> arcs = new ArrayList<>();
    private final Map<String, ChromosomeArc> arcByName = new HashMap<>();
    private final double gapAngle;

    /**
     * @param assembly     the genome whose chromosomes are placed around the circle
     * @param centerX      circle center x in pixels
     * @param centerY      circle center y in pixels
     * @param gapFraction  fraction of the full circle (0..1) reserved, in total,
     *                     for the gaps between chromosomes
     */
    public GenomeArcLayout(Assembly assembly, double centerX, double centerY, double gapFraction) {
        this.centerX = centerX;
        this.centerY = centerY;

        List<Chromosome> chromosomes = assembly.getChromosomes();
        int n = chromosomes.size();
        long totalBp = assembly.totalBp();

        double totalGap = (n > 0) ? gapFraction * 2.0 * Math.PI : 0.0;
        this.gapAngle = (n > 0) ? totalGap / n : 0.0;
        double available = 2.0 * Math.PI - totalGap;

        double angle = START_ANGLE;
        for (Chromosome c : chromosomes) {
            double span = (totalBp > 0) ? available * (c.getBpLength() / (double) totalBp) : 0.0;
            ChromosomeArc arc = new ChromosomeArc(c, angle, angle + span);
            arcs.add(arc);
            arcByName.put(c.getName(), arc);
            angle += span + gapAngle;
        }
    }

    public List<ChromosomeArc> getArcs() {
        return arcs;
    }

    public double getCenterX() {
        return centerX;
    }

    public double getCenterY() {
        return centerY;
    }

    /** Uniform gap angle between adjacent chromosomes, in radians. */
    public double getGapAngle() {
        return gapAngle;
    }

    /**
     * Angle for a base position on a chromosome. Returns {@code NaN} if the
     * chromosome is unknown.
     */
    public double bpToAngle(String refName, long bp) {
        ChromosomeArc arc = arcByName.get(refName);
        if (arc == null) {
            return Double.NaN;
        }
        long len = arc.chromosome.getBpLength();
        double frac = (len > 0) ? Math.max(0.0, Math.min(1.0, bp / (double) len)) : 0.0;
        return arc.startAngle + frac * arc.span();
    }

    /** True if the layout has an arc for the given chromosome name. */
    public boolean hasChromosome(String refName) {
        return arcByName.containsKey(refName);
    }

    /** Point at the given angle and radius, relative to the circle center. */
    public Point2D.Double pointAt(double angle, double radius) {
        return new Point2D.Double(
                centerX + radius * Math.cos(angle),
                centerY + radius * Math.sin(angle));
    }
}
