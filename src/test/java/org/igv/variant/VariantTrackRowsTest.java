package org.igv.variant;

import org.igv.AbstractHeadlessTest;
import org.igv.track.RenderContext;
import org.igv.track.Track;
import org.igv.track.TrackLoader;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Overlapping variants are stacked in rows in EXPANDED and SQUISHED modes, and drawn on a single row in COLLAPSED
 * mode.  The variant found under the mouse -- for tooltips and "show details on click" -- must be the one drawn there.
 */
public class VariantTrackRowsTest extends AbstractHeadlessTest {

    private static final int WIDTH = 1000;

    private File vcf;

    @Before
    public void writeVcf() throws Exception {
        vcf = new File(TestUtils.TMP_OUTPUT_DIR, "overlapping_variants.vcf");
        try (PrintWriter writer = new PrintWriter(vcf)) {
            writer.println("##fileformat=VCFv4.2");
            writer.println("##contig=<ID=chr1,length=8033585>");
            writer.println("##contig=<ID=chr2,length=8033585>");
            writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
            writer.println("chr1\t1000\tdel\tAAAAAAAAAAA\tA\t100\t.\t.");     // Deletes 1001-1010
            writer.println("chr1\t1005\tsnp\tA\tC\t100\t.\t.");               // Overlaps the deletion
            writer.println("chr2\t5000\tsolo\tA\tC\t100\t.\t.");              // Alone, on another chromosome
        }
    }

    @Test
    public void testExpanded() {
        assertVariantUnderMouseIsDrawn(Track.DisplayMode.EXPANDED, 2, 25);
    }

    @Test
    public void testSquished() {
        assertVariantUnderMouseIsDrawn(Track.DisplayMode.SQUISHED, 2, 6);
    }

    @Test
    public void testCollapsed() {
        ReferenceFrame frame = frame();
        VariantTrack track = loadTrack(Track.DisplayMode.COLLAPSED, frame);
        assertEquals("Collapsed draws all variants on one row", 1, track.getNumberOfFeatureLevels());

        BufferedImage image = render(track, frame);
        for (Variant variant : List.of(variant(track, "snp"), variant(track, "del"))) {
            int x = pixel(position(variant), frame);
            for (int y = 0; y < image.getHeight(); y++) {
                if ((image.getRGB(x, y) >>> 24) != 0) {
                    assertTrue("Nothing drawn below the single row, y=" + y, y < 25);
                    assertNotNull("A variant is found where one is drawn, y=" + y,
                            track.getFeatureClosest(position(variant), y, frame, 10 * frame.getScale()));
                }
            }
        }
    }

    /**
     * The deletion covers the SNP's position, so the SNP's column shows the deletion's bar in its row as well as the
     * SNP's bar in its own.  The rows the deletion is drawn in are taken from a column only the deletion covers.
     */
    /**
     * In a multi-panel view (gene list) each panel is packed separately.  A panel with fewer rows than the tallest is
     * still searched without error when the mouse is below its rows.
     */
    @Test
    public void testPanelWithFewerRows() {
        VariantTrack track = loadTrack(Track.DisplayMode.EXPANDED, frame());
        ReferenceFrame soloFrame = new ReferenceFrame("solo");
        soloFrame.setBounds(0, WIDTH);
        soloFrame.jumpTo("chr2", 4990, 5020);
        track.load(soloFrame);
        assertEquals("The tallest panel has two rows", 2, track.getNumberOfFeatureLevels());

        // y = 30 is row 1, which the chr2 panel doesn't have
        Variant found = track.getFeatureClosest(4999.5, 30, soloFrame, 10 * soloFrame.getScale());
        assertNotNull(found);
        assertEquals("solo", found.getID());
    }

    /**
     * A collapsed track doesn't draw genotypes, and "Show Genotypes" says so.  Turning it on expands the track.
     */
    @Test
    public void testShowGenotypesWhenCollapsed() throws Exception {
        File withSamples = new File(TestUtils.TMP_OUTPUT_DIR, "variants_with_samples.vcf");
        try (PrintWriter writer = new PrintWriter(withSamples)) {
            writer.println("##fileformat=VCFv4.2");
            writer.println("##contig=<ID=chr1,length=8033585>");
            writer.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1");
            writer.println("chr1\t1005\tsnp\tA\tC\t100\t.\t.\tGT\t0/1");
        }
        VariantTrack track = (VariantTrack) new TrackLoader().load(new ResourceLocator(withSamples.getAbsolutePath()), genome).get(0);
        assertTrue(track.areGenotypesShown());

        track.setDisplayMode(Track.DisplayMode.COLLAPSED);
        assertFalse(track.areGenotypesShown());

        track.setGenotypesShown(true);
        assertEquals(Track.DisplayMode.EXPANDED, track.getDisplayMode());
        assertTrue(track.areGenotypesShown());

        track.setGenotypesShown(false);
        assertEquals(Track.DisplayMode.EXPANDED, track.getDisplayMode());
        assertFalse(track.areGenotypesShown());
    }

    private void assertVariantUnderMouseIsDrawn(Track.DisplayMode mode, int rows, int rowHeight) {

        ReferenceFrame frame = frame();
        VariantTrack track = loadTrack(mode, frame);
        assertEquals(rows, track.getNumberOfFeatureLevels());
        assertEquals(rowHeight, track.getVariantBandHeight());

        BufferedImage image = render(track, frame);
        Variant snp = variant(track, "snp");
        Variant del = variant(track, "del");

        double deletionOnly = del.getStart() + 1.5;
        List<Integer> deletionRows = drawnRows(image, pixel(deletionOnly, frame));
        List<Integer> snpColumn = drawnRows(image, pixel(position(snp), frame));
        assertTrue("The deletion is drawn", !deletionRows.isEmpty());
        assertTrue("The SNP is drawn below the deletion", snpColumn.stream().anyMatch(y -> !deletionRows.contains(y)));

        List<String> misses = new ArrayList<>();
        for (int y : deletionRows) {
            checkFound(track, frame, deletionOnly, y, "del", misses);
        }
        for (int y : snpColumn) {
            checkFound(track, frame, position(snp), y, deletionRows.contains(y) ? "del" : "snp", misses);
        }
        assertTrue(misses.size() + " misses: " + misses, misses.isEmpty());
    }

    private static void checkFound(VariantTrack track, ReferenceFrame frame, double position, int y, String expected,
                                   List<String> misses) {
        Variant found = track.getFeatureClosest(position, y, frame, 10 * frame.getScale());
        if (found == null || !found.getID().equals(expected)) {
            misses.add(expected + " drawn at y=" + y + ", mouse finds " + (found == null ? "nothing" : found.getID()));
        }
    }

    private static List<Integer> drawnRows(BufferedImage image, int x) {
        List<Integer> ys = new ArrayList<>();
        for (int y = 0; y < image.getHeight(); y++) {
            if ((image.getRGB(x, y) >>> 24) != 0) {
                ys.add(y);
            }
        }
        return ys;
    }

    private static ReferenceFrame frame() {
        ReferenceFrame frame = new ReferenceFrame("test");
        frame.setBounds(0, WIDTH);
        frame.jumpTo("chr1", 990, 1020);
        return frame;
    }

    private VariantTrack loadTrack(Track.DisplayMode mode, ReferenceFrame frame) {
        VariantTrack track = (VariantTrack) new TrackLoader().load(new ResourceLocator(vcf.getAbsolutePath()), genome).get(0);
        track.setDisplayMode(mode);
        track.load(frame);
        return track;
    }

    private static BufferedImage render(VariantTrack track, ReferenceFrame frame) {
        int height = track.getHeight();
        BufferedImage image = new BufferedImage(WIDTH, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D graphics = image.createGraphics();
        Rectangle rect = new Rectangle(0, 0, WIDTH, height);
        track.render(new RenderContext(null, graphics, frame, rect, rect, rect));
        return image;
    }

    private static Variant variant(VariantTrack track, String id) {
        return track.getFeatures("chr1", 0, 2000).stream().map(f -> (Variant) f)
                .filter(v -> id.equals(v.getID())).findFirst().orElseThrow();
    }

    private static double position(Variant variant) {
        return (variant.getStart() + variant.getEnd()) / 2.0;
    }

    private static int pixel(double position, ReferenceFrame frame) {
        return (int) ((position - frame.getOrigin()) / frame.getScale());
    }
}
