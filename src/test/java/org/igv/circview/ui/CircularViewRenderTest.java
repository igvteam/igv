package org.igv.circview.ui;

import org.igv.circview.model.Assembly;
import org.igv.circview.model.Chord;
import org.igv.circview.model.Chromosome;
import org.igv.circview.model.Mate;
import org.igv.circview.render.GenomeArcLayout;
import org.junit.Test;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Headless smoke test: paint into an image, then confirm chord shapes were
 * recorded and that a point on a chord anchor hit-tests to that chord.
 */
public class CircularViewRenderTest {

    private static final int SIZE = 700;

    private static Assembly assembly() {
        return new Assembly("test", "test", List.of(
                new Chromosome("1", 1000, Color.RED),
                new Chromosome("2", 1000, Color.GREEN),
                new Chromosome("3", 1000, Color.BLUE)));
    }

    private CircularView paintedView() {
        CircularViewConfig config = new CircularViewConfig();
        config.width = SIZE;
        config.height = SIZE;
        CircularView view = new CircularView(config);
        view.setSize(SIZE, SIZE);
        view.setAssembly(assembly());
        view.addChords(List.of(
                new Chord("c1", "1", 400, 600, new Mate("2", 400, 600), null)),
                "set", Color.BLUE);

        BufferedImage img = new BufferedImage(SIZE, SIZE, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = img.createGraphics();
        view.paint(g);
        g.dispose();
        return view;
    }

    @Test
    public void paintRecordsChordShapes() {
        CircularView view = paintedView();
        assertEquals(1, view.renderedChordCount());
    }

    @Test
    public void hitTestFindsChordAtAnchor() {
        CircularView view = paintedView();

        // Recompute the geometry the same way the view does to find an on-chord point.
        double outer = SIZE / 2.0 - new CircularViewConfig().margin;
        double ringThickness = outer * new CircularViewConfig().ringThicknessFraction;
        double chordRadius = outer - ringThickness;
        GenomeArcLayout layout = new GenomeArcLayout(assembly(), SIZE / 2.0, SIZE / 2.0,
                new CircularViewConfig().gapFraction);
        // Midpoint of region 1 (bp 500), pulled slightly inward so it lies inside the ribbon.
        double mid = layout.bpToAngle("1", 500);
        Point2D p = layout.pointAt(mid, chordRadius - 2);

        Chord hit = view.chordAt((int) Math.round(p.getX()), (int) Math.round(p.getY()));
        assertNotNull("expected a chord under the anchor point", hit);
        assertEquals("c1", hit.getUniqueId());
    }

    @Test
    public void clickOutsideHitsNothing() {
        CircularView view = paintedView();
        assertNull(view.chordAt(0, 0));
    }

    @Test
    public void tooltipShowsFeatureToStringOverChord() {
        CircularView view = paintedView();

        double outer = SIZE / 2.0 - new CircularViewConfig().margin;
        double chordRadius = outer - outer * new CircularViewConfig().ringThicknessFraction;
        GenomeArcLayout layout = new GenomeArcLayout(assembly(), SIZE / 2.0, SIZE / 2.0,
                new CircularViewConfig().gapFraction);
        Point2D p = layout.pointAt(layout.bpToAngle("1", 500), chordRadius - 2);

        Chord feature = view.chordAt((int) Math.round(p.getX()), (int) Math.round(p.getY()));
        assertNotNull(feature);

        MouseEvent over = new MouseEvent(view, MouseEvent.MOUSE_MOVED, 0L, 0,
                (int) Math.round(p.getX()), (int) Math.round(p.getY()), 0, false);
        assertEquals(feature.tooltipString(), view.getToolTipText(over));
    }

    @Test
    public void tooltipNullAwayFromChords() {
        CircularView view = paintedView();
        MouseEvent off = new MouseEvent(view, MouseEvent.MOUSE_MOVED, 0L, 0, 0, 0, 0, false);
        assertNull(view.getToolTipText(off));
    }

    @Test
    public void hoveredChordIsDrawnBlack() {
        CircularViewConfig config = new CircularViewConfig();
        config.width = SIZE;
        config.height = SIZE;
        CircularView view = new CircularView(config);
        view.setSize(SIZE, SIZE);
        view.setAssembly(assembly());
        Chord c = new Chord("c1", "1", 400, 600, new Mate("2", 400, 600), null);
        view.addChords(List.of(c), "set", Color.BLUE);

        double outer = SIZE / 2.0 - config.margin;
        double chordRadius = outer - outer * config.ringThicknessFraction;
        GenomeArcLayout layout = new GenomeArcLayout(assembly(), SIZE / 2.0, SIZE / 2.0,
                config.gapFraction);
        // A point on the chord just inside the ring (below the chromosome band).
        Point2D p = layout.pointAt(layout.bpToAngle("1", 500), chordRadius - 6);
        int x = (int) Math.round(p.getX());
        int y = (int) Math.round(p.getY());

        // Not hovered: the pixel shows the chord's (blue) color.
        int normal = paintAndSample(view, x, y);
        assertTrue("expected a blue-ish chord pixel, got #" + Integer.toHexString(normal),
                blue(normal) > 120 && red(normal) < 120);

        // Hovered: the same pixel is drawn (near) black.
        view.setHoveredChordForTest(c);
        int hovered = paintAndSample(view, x, y);
        assertTrue("expected a dark hovered pixel, got #" + Integer.toHexString(hovered),
                red(hovered) < 60 && green(hovered) < 60 && blue(hovered) < 60);
    }

    private static int paintAndSample(CircularView view, int x, int y) {
        BufferedImage img = new BufferedImage(SIZE, SIZE, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = img.createGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, SIZE, SIZE);
        view.paint(g);
        g.dispose();
        return img.getRGB(x, y);
    }

    private static int red(int rgb) {
        return (rgb >> 16) & 0xff;
    }

    private static int green(int rgb) {
        return (rgb >> 8) & 0xff;
    }

    private static int blue(int rgb) {
        return rgb & 0xff;
    }

    @Test
    public void overlappingChordsResolveToNearestCenterline() {
        // Two chords share a chr1 anchor but run to chr2 vs chr3. A point at the
        // chr2 anchor must resolve to the chr1->chr2 chord, and the chr3 anchor
        // to the chr1->chr3 chord, even though both pass through the center.
        CircularViewConfig config = new CircularViewConfig();
        config.width = SIZE;
        config.height = SIZE;
        CircularView view = new CircularView(config);
        view.setSize(SIZE, SIZE);
        view.setAssembly(assembly());
        view.addChords(List.of(
                new Chord("to2", "1", 490, 510, new Mate("2", 490, 510), null)),
                "A", Color.BLUE);
        view.addChords(List.of(
                new Chord("to3", "1", 490, 510, new Mate("3", 490, 510), null)),
                "B", Color.RED);
        BufferedImage img = new BufferedImage(SIZE, SIZE, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = img.createGraphics();
        view.paint(g);
        g.dispose();

        double outer = SIZE / 2.0 - config.margin;
        double chordRadius = outer - outer * config.ringThicknessFraction;
        GenomeArcLayout layout = new GenomeArcLayout(assembly(), SIZE / 2.0, SIZE / 2.0,
                config.gapFraction);

        Point2D atChr2 = layout.pointAt(layout.bpToAngle("2", 500), chordRadius);
        Chord hit2 = view.chordAt((int) Math.round(atChr2.getX()), (int) Math.round(atChr2.getY()));
        assertNotNull(hit2);
        assertEquals("to2", hit2.getUniqueId());

        Point2D atChr3 = layout.pointAt(layout.bpToAngle("3", 500), chordRadius);
        Chord hit3 = view.chordAt((int) Math.round(atChr3.getX()), (int) Math.round(atChr3.getY()));
        assertNotNull(hit3);
        assertEquals("to3", hit3.getUniqueId());
    }
}
