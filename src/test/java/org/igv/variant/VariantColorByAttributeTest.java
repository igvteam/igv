package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.igv.AbstractHeadlessTest;
import org.igv.track.RenderContext;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.Before;
import org.junit.Test;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests for coloring the variant band by a VCF INFO attribute (issue #1657).
 */
public class VariantColorByAttributeTest extends AbstractHeadlessTest {

    private VariantTrack track;
    private List<Feature> variants;

    @Before
    public void loadTrack() throws Exception {
        String filePath = TestUtils.DATA_DIR + "vcf/clinvar_info.vcf";
        TestUtils.createIndex(filePath);
        track = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(filePath), genome).get(0);
        variants = track.getFeatures("chr1", 0, 1000);
        assertEquals(9, variants.size());
    }

    /**
     * Only categorical INFO fields are offered, and they are sorted by ID.  AF is Float, hence excluded.
     */
    @Test
    public void testColorableInfoFields() {
        List<String> ids = track.getColorableInfoFields().stream()
                .map(VCFInfoHeaderLine::getID)
                .collect(Collectors.toList());
        assertEquals(List.of("ALLELEID", "CLNREVSTAT", "CLNSIG", "DB", "SVTYPE"), ids);
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
     * A Float attribute has no default treatment -- it is only offered if a scheme gives it a range.
     */
    @Test
    public void testFloatAttributeNotOffered() {
        List<String> ids = track.getColorableInfoFields().stream()
                .map(VCFInfoHeaderLine::getID)
                .collect(Collectors.toList());
        assertFalse("AF is a Float attribute with no scheme", ids.contains("AF"));
    }

    /**
     * Coloring is categorical, so an attribute with unbounded values must not get a color per variant.  This is
     * reachable even though the menu does not offer Float attributes -- a session or a batch command can name
     * any attribute, and an Integer or String attribute can be just as unbounded.
     */
    @Test
    public void testColorLimit() {
        int limit = VariantTrack.getMaxAttributeColors();
        track.setColorByAttribute("SCORE");
        assertFalse(track.isAttributeColorLimitReached("SCORE"));

        for (int i = 0; i < limit + 20; i++) {
            track.getAttributeColor("SCORE", "value-" + i);
        }

        assertTrue(track.isAttributeColorLimitReached("SCORE"));
        assertEquals(limit, track.getAttributeColorTable("SCORE").getColorMap().size());
        assertEquals(Color.gray, track.getAttributeColor("SCORE", "value-" + (limit + 100)));
    }

    private static double distance(Color c1, Color c2) {
        int dr = c1.getRed() - c2.getRed();
        int dg = c1.getGreen() - c2.getGreen();
        int db = c1.getBlue() - c2.getBlue();
        return Math.sqrt(dr * dr + dg * dg + db * db);
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
