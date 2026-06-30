package org.igv.circview.ui;

import org.igv.circview.model.Assembly;
import org.igv.circview.model.ChordCollection;
import org.igv.circview.model.Chord;
import org.igv.circview.model.ChordSet;
import org.igv.circview.model.ChordSetManager;
import org.igv.circview.model.Mate;
import org.igv.circview.render.GenomeArcLayout;

import javax.swing.JPanel;
import javax.swing.ToolTipManager;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

/**
 * A Swing component that renders a genome as a circle and draws chords between
 * interacting regions (a "Circos" view).
 *
 * <p>This is a clean Java2D reimplementation of the view that circularView.js
 * delegates to JBrowse. It renders the chromosome ring and chords, hit-tests
 * clicks (firing {@link CircularViewConfig#onChordClick}), and shows a hover
 * tooltip for the chord under the cursor. The surrounding toolbar and control
 * panel live in {@link CircularViewPanel}.
 */
public class CircularView extends JPanel {

    /**
     * A chord's rendered geometry, retained for hit-testing. Stores both the
     * filled ribbon (for "inside the ribbon" hits on wide chords) and the
     * chord's centerline — the quadratic Bézier from region A's midpoint,
     * through the circle center, to region B's midpoint — used to pick the
     * nearest chord among overlapping ones.
     */
    private static final class RenderedChord {
        final Path2D shape;
        final Rectangle2D bounds;
        final Chord feature;
        // Centerline (spine): quadratic Bézier (sx0,sy0) -ctrl(scx,scy)-> (sx2,sy2).
        final double sx0, sy0, scx, scy, sx2, sy2;

        RenderedChord(Path2D shape, Chord feature,
                      double sx0, double sy0, double scx, double scy, double sx2, double sy2) {
            this.shape = shape;
            this.bounds = shape.getBounds2D();
            this.feature = feature;
            this.sx0 = sx0;
            this.sy0 = sy0;
            this.scx = scx;
            this.scy = scy;
            this.sx2 = sx2;
            this.sy2 = sy2;
        }

        /** Squared distance from (px,py) to the centerline, sampled along the curve. */
        double centerlineDistanceSq(double px, double py) {
            final int segments = 16;
            double best = Double.MAX_VALUE;
            double prevX = sx0, prevY = sy0;
            for (int i = 1; i <= segments; i++) {
                double t = i / (double) segments;
                double mt = 1 - t;
                double bx = mt * mt * sx0 + 2 * mt * t * scx + t * t * sx2;
                double by = mt * mt * sy0 + 2 * mt * t * scy + t * t * sy2;
                double d = segmentDistanceSq(px, py, prevX, prevY, bx, by);
                if (d < best) {
                    best = d;
                }
                prevX = bx;
                prevY = by;
            }
            return best;
        }
    }

    /** Slightly heavier stroke for the hovered chord so thin arcs stand out. */
    private static final BasicStroke HOVER_STROKE = new BasicStroke(1.75f);

    private final CircularViewConfig config;
    private final ChordSetManager chordManager = new ChordSetManager();
    private Assembly assembly;
    private boolean groupByTrack = false;

    private final List<RenderedChord> renderedChords = new ArrayList<>();

    /** The chord currently under the cursor; drawn highlighted. */
    private Chord hoveredChord;

    /** Notified when the set of collections changes (added/cleared/regrouped). */
    private final List<Runnable> structureListeners = new ArrayList<>();

    public CircularView(CircularViewConfig config) {
        this.config = (config != null) ? config : new CircularViewConfig();
        setBackground(Color.WHITE);
        setPreferredSize(new Dimension(this.config.width, this.config.height));
        MouseAdapter mouse = new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                handleClick(e.getX(), e.getY());
            }

            @Override
            public void mouseMoved(MouseEvent e) {
                setHoveredChord(chordAt(e.getX(), e.getY()));
            }

            @Override
            public void mouseExited(MouseEvent e) {
                setHoveredChord(null);
            }
        };
        addMouseListener(mouse);
        addMouseMotionListener(mouse);
        // Enable per-position tooltips (see getToolTipText below).
        ToolTipManager.sharedInstance().registerComponent(this);
    }

    /** Update the highlighted chord, repainting only when it actually changes. */
    private void setHoveredChord(Chord feature) {
        if (feature != hoveredChord) {
            hoveredChord = feature;
            repaint();
        }
    }

    /** Tooltip text for the chord under the cursor, or null when over none. */
    @Override
    public String getToolTipText(MouseEvent event) {
        Chord feature = chordAt(event.getX(), event.getY());
        return (feature != null) ? feature.tooltipString() : null;
    }

    // ---- Public API (subset of the JS CircularView) -------------------------

    /** Reset the view with a new genome. */
    public void setAssembly(Assembly assembly) {
        this.assembly = assembly;
        chordManager.clearChords();
        fireStructureChanged();
        repaint();
    }

    /**
     * Append (or replace by name) a set of chords.
     *
     * @param chords     the chord features
     * @param name       chord-set name; the track name is the substring before
     *                   the first space, matching addChords() in circularView.js
     * @param color      chord-set color (used when a feature has no color)
     * @param trackColor track color; defaults to {@code color} when null
     */
    public void addChords(List<Chord> chords, String name, Color color, Color trackColor) {
        String setName = (name != null) ? name : "*";
        String trackName = setName.split(" ")[0];
        Color c = (color != null) ? color : Color.BLACK;
        Color tc = (trackColor != null) ? trackColor : c;
        chordManager.addChordSet(new ChordSet(setName, trackName, chords, c, tc));
        fireStructureChanged();
        repaint();
    }

    public void addChords(List<Chord> chords, String name, Color color) {
        addChords(chords, name, color, color);
    }

    public void clearChords() {
        chordManager.clearChords();
        fireStructureChanged();
        repaint();
    }

    public boolean isGroupByTrack() {
        return groupByTrack;
    }

    public void setGroupByTrack(boolean groupByTrack) {
        this.groupByTrack = groupByTrack;
        fireStructureChanged();
        repaint();
    }

    /**
     * The collections currently driving the view: the flat chord sets, or the
     * tracks when grouping by track. Mirrors getChordSet() in circularView.js,
     * which switches on the groupByTrack flag.
     */
    public List<? extends ChordCollection> getActiveCollections() {
        return groupByTrack ? chordManager.getTracks() : chordManager.getChordSets();
    }

    private ChordCollection activeByName(String name) {
        for (ChordCollection c : getActiveCollections()) {
            if (c.getName().equals(name)) {
                return c;
            }
        }
        return null;
    }

    /** Set the color of the active collection (chord set or track) with this name. */
    public void setColor(String name, Color color) {
        ChordCollection c = activeByName(name);
        if (c != null) {
            c.setColor(color);
            repaint();
        }
    }

    public void hideChordSet(String name) {
        setVisible(name, false);
    }

    public void showChordSet(String name) {
        setVisible(name, true);
    }

    public void setVisible(String name, boolean visible) {
        ChordCollection c = activeByName(name);
        if (c != null) {
            c.setVisible(visible);
            repaint();
        }
    }

    public ChordSetManager getChordManager() {
        return chordManager;
    }

    public CircularViewConfig getConfig() {
        return config;
    }

    /** Register a listener fired when collections are added, cleared, or regrouped. */
    public void addStructureListener(Runnable listener) {
        structureListeners.add(listener);
    }

    private void fireStructureChanged() {
        for (Runnable r : structureListeners) {
            r.run();
        }
    }

    // ---- Rendering ----------------------------------------------------------

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        renderedChords.clear();
        if (assembly == null || assembly.getChromosomes().isEmpty()) {
            return;
        }

        Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

            int w = getWidth();
            int h = getHeight();
            double size = Math.min(w, h);
            double cx = w / 2.0;
            double cy = h / 2.0;

            double outerRadius = size / 2.0 - config.margin;
            if (outerRadius <= 0) {
                return;
            }
            double ringThickness = outerRadius * config.ringThicknessFraction;
            double ringInner = outerRadius - ringThickness;
            double chordRadius = ringInner;

            GenomeArcLayout layout = new GenomeArcLayout(assembly, cx, cy, config.gapFraction);

            drawChords(g2, layout, chordRadius, cx, cy);
            drawRing(g2, layout, ringInner, outerRadius);
            if (config.showLabels) {
                drawLabels(g2, layout, outerRadius);
            }
        } finally {
            g2.dispose();
        }
    }

    private void drawRing(Graphics2D g2, GenomeArcLayout layout, double ringInner, double ringOuter) {
        for (GenomeArcLayout.ChromosomeArc arc : layout.getArcs()) {
            Path2D band = ringBand(layout, arc.startAngle, arc.endAngle, ringInner, ringOuter);
            g2.setColor(arc.chromosome.getColor());
            g2.fill(band);
        }
    }

    private void drawLabels(Graphics2D g2, GenomeArcLayout layout, double outerRadius) {
        g2.setColor(Color.DARK_GRAY);
        g2.setFont(getFont().deriveFont(Font.BOLD, 11f));
        FontMetrics fm = g2.getFontMetrics();
        double labelRadius = outerRadius + 14;
        for (GenomeArcLayout.ChromosomeArc arc : layout.getArcs()) {
            double mid = (arc.startAngle + arc.endAngle) / 2.0;
            Point2D p = layout.pointAt(mid, labelRadius);
            String name = arc.chromosome.getName();
            int tw = fm.stringWidth(name);
            g2.drawString(name, (float) (p.getX() - tw / 2.0),
                    (float) (p.getY() + fm.getAscent() / 2.0 - 1));
        }
    }

    private void drawChords(Graphics2D g2, GenomeArcLayout layout, double chordRadius, double cx, double cy) {
        g2.setStroke(new BasicStroke(0.75f));
        RenderedChord hovered = null;
        for (ChordCollection collection : getActiveCollections()) {
            if (!collection.isVisible()) {
                continue;
            }
            for (Chord chord : collection.getChords()) {
                Mate mate = chord.getMate();
                if (mate == null
                        || !layout.hasChromosome(chord.getRefName())
                        || !layout.hasChromosome(mate.getRefName())) {
                    continue;
                }
                double a0 = layout.bpToAngle(chord.getRefName(), chord.getStart());
                double a1 = layout.bpToAngle(chord.getRefName(), chord.getEnd());
                double b0 = layout.bpToAngle(mate.getRefName(), mate.getStart());
                double b1 = layout.bpToAngle(mate.getRefName(), mate.getEnd());

                Path2D ribbon = chordRibbon(layout, a0, a1, b0, b1, chordRadius, cx, cy);
                Color color = (chord.getColor() != null) ? chord.getColor() : collection.getColor();
                // Fill and outline both use the chord's alpha so the transparency
                // slider is honored. (An opaque outline would mask the fill alpha,
                // which matters because most chords are very thin ribbons.)
                g2.setColor(color);
                g2.fill(ribbon);
                g2.draw(ribbon);

                Point2D aMid = layout.pointAt((a0 + a1) / 2.0, chordRadius);
                Point2D bMid = layout.pointAt((b0 + b1) / 2.0, chordRadius);
                RenderedChord rc = new RenderedChord(ribbon, chord,
                        aMid.getX(), aMid.getY(), cx, cy, bMid.getX(), bMid.getY());
                renderedChords.add(rc);
                if (chord == hoveredChord) {
                    hovered = rc;
                }
            }
        }
        // Redraw the hovered chord opaque-black and on top, so it is unambiguous
        // among overlapping arcs before the user clicks.
        if (hovered != null) {
            g2.setColor(Color.BLACK);
            g2.setStroke(HOVER_STROKE);
            g2.fill(hovered.shape);
            g2.draw(hovered.shape);
        }
    }

    /**
     * A ribbon connecting region [start,end] on one chromosome to [start,end] on
     * its mate, with both cross-edges bowing toward the circle center.
     */
    private Path2D chordRibbon(GenomeArcLayout layout, double a0, double a1, double b0, double b1,
                               double r, double cx, double cy) {
        Point2D pA0 = layout.pointAt(a0, r);
        Point2D pB1 = layout.pointAt(b1, r);

        Path2D path = new Path2D.Double();
        path.moveTo(pA0.getX(), pA0.getY());
        appendRingArc(path, layout, a0, a1, r);                 // along ring, region A
        path.quadTo(cx, cy, pB1.getX(), pB1.getY());            // cross to region B (bow to center)
        appendRingArc(path, layout, b1, b0, r);                 // along ring, region B
        path.quadTo(cx, cy, pA0.getX(), pA0.getY());            // cross back to region A
        path.closePath();
        return path;
    }

    /** Append points along the ring from angle a to angle b (inclusive of b). */
    private void appendRingArc(Path2D path, GenomeArcLayout layout, double a, double b, double r) {
        double sweep = b - a;
        int steps = Math.max(1, (int) Math.ceil(Math.abs(sweep) / 0.05));
        for (int i = 1; i <= steps; i++) {
            double t = a + sweep * (i / (double) steps);
            Point2D p = layout.pointAt(t, r);
            path.lineTo(p.getX(), p.getY());
        }
    }

    private Path2D ringBand(GenomeArcLayout layout, double a, double b, double rInner, double rOuter) {
        Path2D path = new Path2D.Double();
        Point2D start = layout.pointAt(a, rOuter);
        path.moveTo(start.getX(), start.getY());
        appendRingArc(path, layout, a, b, rOuter);   // outer edge a -> b
        Point2D inner = layout.pointAt(b, rInner);
        path.lineTo(inner.getX(), inner.getY());
        appendRingArc(path, layout, b, a, rInner);   // inner edge b -> a
        path.closePath();
        return path;
    }

    // ---- Hit-testing --------------------------------------------------------

    private void handleClick(int x, int y) {
        Chord feature = chordAt(x, y);
        if (feature != null && config.onChordClick != null) {
            config.onChordClick.onChordClick(feature);
        }
    }

    /** Maximum distance, in pixels, from a chord's centerline to count as a hit. */
    private static final int HIT_TOLERANCE = 4;

    /**
     * The chord under the point, or null. Among overlapping chords the one whose
     * centerline passes nearest the point is chosen; a point inside a (wide)
     * ribbon counts as distance zero. The match must lie within
     * {@link #HIT_TOLERANCE} pixels. On exact ties the top-most (last drawn)
     * chord wins. Linear in the number of chords, with a bounding-box prune.
     *
     * <p>Valid only after a paint has populated the rendered shapes.
     */
    Chord chordAt(int x, int y) {
        final double tol = HIT_TOLERANCE;
        final double tolSq = tol * tol;
        double bestSq = Double.MAX_VALUE;
        Chord best = null;
        for (RenderedChord rc : renderedChords) {
            // Cheap reject: a point outside the ribbon's bounds (plus tolerance)
            // can't be within tolerance of its centerline either.
            Rectangle2D b = rc.bounds;
            if (x < b.getMinX() - tol || x > b.getMaxX() + tol
                    || y < b.getMinY() - tol || y > b.getMaxY() + tol) {
                continue;
            }
            double dSq = rc.shape.contains(x, y) ? 0.0 : rc.centerlineDistanceSq(x, y);
            if (dSq <= bestSq) { // <= so a later (top-most) chord wins ties
                bestSq = dSq;
                best = rc.feature;
            }
        }
        return (bestSq <= tolSq) ? best : null;
    }

    /** Squared distance from point (px,py) to segment (x1,y1)-(x2,y2). */
    private static double segmentDistanceSq(double px, double py,
                                            double x1, double y1, double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        double lenSq = dx * dx + dy * dy;
        double t = (lenSq == 0.0) ? 0.0 : ((px - x1) * dx + (py - y1) * dy) / lenSq;
        t = Math.max(0.0, Math.min(1.0, t));
        double cx = x1 + t * dx;
        double cy = y1 + t * dy;
        double ex = px - cx;
        double ey = py - cy;
        return ex * ex + ey * ey;
    }

    /** Number of chord shapes recorded by the last paint (test hook). */
    int renderedChordCount() {
        return renderedChords.size();
    }

    /** Set the highlighted chord without repainting (test hook). */
    void setHoveredChordForTest(Chord feature) {
        this.hoveredChord = feature;
    }
}
