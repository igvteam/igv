package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.igv.AbstractHeadlessTest;
import org.igv.DirectoryManager;
import org.igv.track.RenderContext;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Tests for coloring the variant band by a VCF INFO attribute (issue #1657).
 */
public class VariantColorByAttributeTest extends AbstractHeadlessTest {

    private VariantTrack track;
    private List<Feature> variants;
    private File previousIgvDirectory;

    /**
     * These tests assert exact built-in colors, so they must not see schemes the developer running them has
     * installed, and must not leave cached schemes behind for the next test class.
     */
    @After
    public void restoreIgvDirectory() {
        DirectoryManager.setIgvDirectory(previousIgvDirectory);
        VariantColorSchemes.reset();
    }

    // One @Before -- JUnit does not order them, and the directory has to be set before anything reads a scheme
    @Before
    public void loadTrack() throws Exception {

        previousIgvDirectory = DirectoryManager.getIgvDirectory();
        File igvDirectory = new File(TestUtils.TMP_OUTPUT_DIR, "igv");
        igvDirectory.mkdirs();
        DirectoryManager.setIgvDirectory(igvDirectory);
        VariantColorSchemes.reset();

        String filePath = TestUtils.DATA_DIR + "vcf/clinvar_info.vcf";
        TestUtils.createIndex(filePath);
        track = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(filePath), genome).get(0);
        variants = track.getFeatures("chr1", 0, 1000);
        assertEquals(9, variants.size());
    }

    /**
     * INFO fields are offered sorted by ID.  Numeric fields are included -- selecting one asks for a color scale.
     */
    @Test
    public void testColorableInfoFields() {
        List<String> ids = track.getColorableInfoFields().stream()
                .map(VCFInfoHeaderLine::getID)
                .collect(Collectors.toList());
        assertEquals(List.of("AF", "ALLELEID", "CLNREVSTAT", "CLNSIG", "DB", "RDP", "SVTYPE"), ids);
    }

    @Test
    public void testSetColorByAttribute() {
        track.setColorByAttribute("CLNSIG");
        assertEquals(VariantTrack.ColorMode.ATTRIBUTE, track.getSiteColorMode());
        assertEquals("CLNSIG", track.getColorByAttribute());

        track.setColorByAttribute(null);
        assertEquals(VariantTrack.ColorMode.NONE, track.getSiteColorMode());
        assertEquals(null, track.getColorByAttribute());
    }

    /**
     * ClinVar significance values have predefined colors, benign (blue) through pathogenic (red).
     */
    @Test
    public void testClinicalSignificanceColors() {
        track.setColorByAttribute("CLNSIG");
        assertEquals(new Color(202, 0, 32), colorAt(0));    // Pathogenic
        assertEquals(new Color(244, 109, 67), colorAt(1));  // Likely_pathogenic
        assertEquals(new Color(150, 150, 150), colorAt(2)); // Uncertain_significance
        assertEquals(new Color(146, 197, 222), colorAt(3)); // Likely_benign
        assertEquals(new Color(5, 113, 176), colorAt(4));   // Benign
    }

    /**
     * Values with no predefined color get distinct colors from the palette, and a variant with no value for the
     * attribute is drawn gray.
     */
    @Test
    public void testUnrecognizedAndMissingValues() {
        track.setColorByAttribute("CLNSIG");
        Color drugResponse = colorAt(5);
        Color association = colorAt(6);
        assertNotEquals(drugResponse, association);
        assertEquals(Color.gray, colorAt(7));       // no CLNSIG attribute
        assertEquals(drugResponse, colorAt(5));     // assignments are stable
    }

    @Test
    public void testStructuralVariantTypeColors() {
        track.setColorByAttribute("SVTYPE");
        assertEquals(new Color(255, 33, 1), colorAt(8));   // DEL
        assertEquals(Color.gray, colorAt(0));              // no SVTYPE attribute
    }

    /**
     * Multi-valued attributes are keyed by their values joined with a comma, not by htsjdk's list rendering
     * ("[a, b]").
     */
    @Test
    public void testMultiValuedAttribute() {
        track.setColorByAttribute("CLNREVSTAT");
        colorAt(0);
        assertTrue(track.getAttributeColorTable("CLNREVSTAT").getKeys()
                .contains("criteria_provided,_single_submitter"));
    }

    /**
     * The selected attribute and its color assignments survive a session round trip -- palette colors are
     * assigned in the order values are seen, so they have to be persisted to be reproducible.
     */
    @Test
    public void testSessionRoundTrip() {
        track.setColorByAttribute("CLNSIG");
        Color drugResponse = colorAt(5);

        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertEquals("CLNSIG", json.getString("colorByAttribute"));
        assertEquals("ATTRIBUTE", json.getString("siteColorMode"));

        VariantTrack restored = new VariantTrack();
        restored.unmarshalJSON(json);
        assertEquals("CLNSIG", restored.getColorByAttribute());
        assertEquals(VariantTrack.ColorMode.ATTRIBUTE, restored.getSiteColorMode());
        assertEquals(drugResponse, restored.getAttributeColor((Variant) variants.get(5)));
        assertFalse(json.getString("attributeColorTable").isEmpty());
    }

    /**
     * A value with no scheme entry must not come out looking like one that has one -- Set 1's red is close enough
     * to the ClinVar "Pathogenic" red to be mistaken for it at the size of a variant band.
     */
    @Test
    public void testPaletteColorsAvoidSchemeColors() {
        track.setColorByAttribute("CLNSIG");

        for (int index : new int[]{5, 6}) {          // drug_response, association -- no scheme entry
            Color assigned = colorAt(index);
            for (Color schemeColor : VariantColorSchemes.getColors("CLNSIG")) {
                assertTrue("Assigned color " + assigned + " is too close to scheme color " + schemeColor,
                        distance(assigned, schemeColor) >= 60);
            }
        }
    }

    /**
     * Colors keep their distance once the palette runs out.  The palette has nine entries and the built-in
     * CLNSIG scheme rejects four of them, so from the sixth unknown value on, colors have to be generated --
     * a fallback that took palette[n] regardless would hand back near-duplicates and repeats.
     */
    @Test
    public void testAssignedColorsStayDistinctBeyondThePalette() {

        track.setColorByAttribute("CLNSIG");

        List<Color> assigned = new ArrayList<>();
        for (int i = 0; i < 25; i++) {
            assigned.add(track.getAttributeColor("CLNSIG", "unknown-" + i));
        }

        Collection<Color> schemeColors = VariantColorSchemes.getColors("CLNSIG");
        for (int i = 0; i < assigned.size(); i++) {
            Color color = assigned.get(i);

            // Beyond a couple of dozen colors the RGB cube runs out of room, so the guarantee is the best
            // available rather than the full separation -- but it must never collapse
            double required = i < 15 ? 60 : 55;

            for (Color schemeColor : schemeColors) {
                assertTrue("Value " + i + " got " + color + ", only " + distance(color, schemeColor)
                                + " from scheme color " + schemeColor,
                        distance(color, schemeColor) >= required);
            }
            for (int j = 0; j < i; j++) {
                assertTrue("Value " + i + " repeats the color of value " + j,
                        distance(color, assigned.get(j)) > 0);
            }
        }
    }

    /**
     * The generated sequence must not cycle -- ColorUtilities.randomColor repeats every 215 indices, and a
     * search over a cycling sequence can only find duplicates once the cycle is used up.  Pairwise distinctness
     * of the candidates is what makes the exhaustion policy exact.
     */
    @Test
    public void testGeneratedColorsDoNotRepeat() {
        java.util.Set<Color> seen = new java.util.HashSet<>();
        for (int i = 0; i < 2000; i++) {
            assertTrue("Generated color " + i + " repeats an earlier one", seen.add(VariantTrack.generatedColor(i)));
        }
    }

    /**
     * Well past the 215 colors the old generator could ever produce, values still get colors of their own.
     */
    @Test
    public void testNoRepeatsBeyondTheOldGeneratorsPeriod() {
        track.setColorByAttribute("CLNSIG");
        java.util.Set<Color> assigned = new java.util.HashSet<>();
        for (int i = 0; i < 400; i++) {
            Color color = track.getAttributeColor("CLNSIG", "unknown-" + i);
            assertTrue("Value " + i + " was given a color already in use", assigned.add(color));
        }
    }

    /**
     * The same value keeps its color, however many others have been assigned since.
     */
    @Test
    public void testAssignedColorsAreStable() {
        track.setColorByAttribute("CLNSIG");
        Color first = track.getAttributeColor("CLNSIG", "unknown-0");
        for (int i = 1; i < 20; i++) {
            track.getAttributeColor("CLNSIG", "unknown-" + i);
        }
        assertEquals(first, track.getAttributeColor("CLNSIG", "unknown-0"));
    }

    private static double distance(Color c1, Color c2) {
        int dr = c1.getRed() - c2.getRed();
        int dg = c1.getGreen() - c2.getGreen();
        int db = c1.getBlue() - c2.getBlue();
        return Math.sqrt(dr * dr + dg * dg + db * db);
    }

    /**
     * "." is the VCF missing value marker.  htsjdk passes it through for a string attribute, so it has to be
     * recognized here or it takes a color of its own and shows in the legend as though it were a real value.
     */
    @Test
    public void testMissingValueMarker() {
        track.setColorByAttribute("CLNREVSTAT");

        assertEquals("CLNREVSTAT=. must read as missing", Color.gray, colorAt(1));
        assertEquals("A record with no CLNREVSTAT at all is missing too", Color.gray, colorAt(2));
        assertFalse(track.getAttributeColorTable("CLNREVSTAT").getColorMap().containsKey("."));
    }

    /**
     * Only "." and an empty value mean missing.  An absent attribute arrives as a Java null, so the text
     * "null" is a value the file actually contains, and must be colorable and visible in the legend.
     */
    @Test
    public void testLiteralNullIsAValue() {
        assertEquals("null", VariantTrack.normalizeAttributeValue("null"));
        assertNull(VariantTrack.normalizeAttributeValue(null));
        assertNull(VariantTrack.normalizeAttributeValue("."));
        assertNull(VariantTrack.normalizeAttributeValue(""));
        assertNull(VariantTrack.normalizeAttributeValue("  "));
        assertEquals("a,b", VariantTrack.normalizeAttributeValue("[a, b]"));
    }

    /**
     * The renderer fills the whole variant band with the attribute color -- no allele frequency bar.
     */
    @Test
    public void testRenderSiteBand() {
        track.setColorByAttribute("CLNSIG");

        BufferedImage image = new BufferedImage(20, 25, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        Rectangle rect = new Rectangle(0, 0, 20, 25);
        RenderContext context = new RenderContext(null, graphics, null, rect, rect, rect);

        new VariantRenderer(track).renderSiteBand((Variant) variants.get(0), rect, 0, 20, context);

        // Pathogenic, top margin is 3 pixels
        assertEquals(new Color(202, 0, 32), new Color(image.getRGB(10, 4)));
        assertEquals(new Color(202, 0, 32), new Color(image.getRGB(10, 24)));
    }

    private Color colorAt(int index) {
        return track.getAttributeColor((Variant) variants.get(index));
    }
}
