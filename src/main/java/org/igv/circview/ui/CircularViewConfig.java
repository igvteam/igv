package org.igv.circview.ui;

/**
 * Rendering and behavior options for a {@link CircularView}. Defaults give a
 * reasonable Circos-style look; all fields may be overridden before the view is
 * shown.
 */
public final class CircularViewConfig {

    /** Nominal view size in pixels (square). */
    public int width = 700;
    public int height = 700;

    /** Fraction of the full circle reserved, in total, for gaps between chromosomes. */
    public double gapFraction = 0.015;

    /** Thickness of the chromosome ring as a fraction of the circle radius. */
    public double ringThicknessFraction = 0.03;

    /** Pixel margin between the circle and the panel edge (room for labels). */
    public int margin = 40;

    /** Whether to draw chromosome name labels outside the ring. */
    public boolean showLabels = true;

    /** Invoked when a chord is clicked. Defaults to printing the feature. */
    public ChordClickListener onChordClick = feature -> System.out.println(feature);
}
