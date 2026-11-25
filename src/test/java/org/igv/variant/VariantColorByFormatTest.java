package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.igv.AbstractHeadlessTest;
import org.igv.renderer.ContinuousColorScale;
import org.igv.prefs.PreferencesManager;
import org.igv.track.DataType;
import org.igv.track.RenderContext;
import org.igv.track.TrackLoader;
import org.igv.ui.panel.ReferenceFrame;
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
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Tests for coloring genotypes by a VCF FORMAT field.
 */
public class VariantColorByFormatTest extends AbstractHeadlessTest {

    private static final ContinuousColorScale GQ_SCALE = new ContinuousColorScale(10, 99, Color.white, new Color(200, 0, 0));
    private static final ContinuousColorScale AD_SCALE = new ContinuousColorScale(0, 9, Color.white, new Color(200, 0, 0));

    private VariantTrack track;
    private List<Feature> variants;

    @Before
    public void loadTrack() throws Exception {
        String filePath = TestUtils.DATA_DIR + "vcf/format_colors.vcf";
        TestUtils.createIndex(filePath);
        track = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(filePath), genome).get(0);
        variants = track.getFeatures("chr1", 0, 1000);
        assertEquals(2, variants.size());
    }

    /**
     * FORMAT fields are offered sorted by ID, without GT -- coloring by genotype is its own option.
     */
    @Test
    public void testColorableFormatFields() {
        List<String> ids = track.getColorableFormatFields().stream()
                .map(VCFFormatHeaderLine::getID)
                .collect(Collectors.toList());
        assertEquals(List.of("AD", "CL", "CN", "GQ"), ids);
    }

    /**
     * What the menu and the default scale read from the loaded features: which declared fields have values at
     * all -- a header can declare fields no record carries -- and the values and range of one that does.
     */
    @Test
    public void testLoadedFeatureValues() {
        ReferenceFrame frame = new ReferenceFrame("test");
        frame.setBounds(0, 500);
        frame.jumpTo("chr1", 50, 250);
        track.load(frame);

        assertEquals(Set.of("AD", "CL", "GQ"), track.getFormatFieldsWithValues(List.of("AD", "CL", "CN", "GQ")));
        assertEquals(List.of("high", "low"), List.copyOf(track.getAttributeValues(VariantTrack.formatKey("CL"))));
        assertArrayEquals(new double[]{10, 99}, track.getAttributeRange(VariantTrack.formatKey("GQ")), 0);
        assertNull(track.getAttributeRange(VariantTrack.formatKey("CN")));
    }

    /**
     * A numeric field is colored by the track's scale.  A genotype with no value has no color.
     */
    @Test
    public void testNumericField() {
        track.setColorScale(VariantTrack.formatKey("GQ"), GQ_SCALE);
        track.setColorGenotypesBy("GQ");
        assertEquals(VariantTrack.ColorMode.FORMAT_ATTRIBUTE, track.getGenotypeColorMode());

        assertEquals(GQ_SCALE.getColor(10), colorOf(0, "S1"));
        assertEquals(GQ_SCALE.getColor(99), colorOf(0, "S3"));
        assertNotEquals(colorOf(0, "S1"), colorOf(0, "S3"));
        assertNull(colorOf(0, "S4"));           // GQ missing
    }

    /**
     * A Number=R field lists the reference allele first.  It is skipped, and the alternates reduced to their maximum.
     */
    @Test
    public void testPerAlleleField() {
        track.setColorScale(VariantTrack.formatKey("AD"), AD_SCALE);
        track.setColorGenotypesBy("AD");

        assertEquals(AD_SCALE.getColor(9), colorOf(1, "S1"));    // 3,4,9
        assertEquals(AD_SCALE.getColor(0), colorOf(1, "S3"));    // 15,0,0
    }

    /**
     * A string field gets a color per value.  The same ID as an INFO attribute is colored separately.
     */
    @Test
    public void testCategoricalField() {
        track.setColorGenotypesBy("CL");

        assertEquals(colorOf(0, "S2"), colorOf(0, "S3"));       // high, high
        assertNotEquals(colorOf(0, "S1"), colorOf(0, "S2"));    // low, high
        assertNull(colorOf(0, "S4"));

        Color chosen = new Color(1, 2, 3);
        track.setAttributeColor(VariantTrack.formatKey("CL"), "low", chosen);
        assertEquals(chosen, colorOf(0, "S1"));
        assertNotEquals(chosen, track.getAttributeColor("CL", "low"));
    }

    /**
     * CN defaults to the colors of segmented copy number tracks: neutral at 2 copies, blue below, red above.
     */
    @Test
    public void testCopyNumberDefaultScale() {
        String key = VariantTrack.formatKey("CN");
        assertTrue(VariantTrackMenuHelper.defineScaleIfNeeded(track, key));

        ContinuousColorScale segScale = PreferencesManager.getPreferences().getColorScale(DataType.COPY_NUMBER);
        assertClose(segScale.getMinColor(), track.getAttributeColor(key, "0"));
        assertClose(segScale.getMidColor(), track.getAttributeColor(key, "2"));
        assertClose(segScale.getMaxColor(), track.getAttributeColor(key, "6"));
        assertClose(segScale.getMaxColor(), track.getAttributeColor(key, "9"));
        assertNotEquals(track.getAttributeColor(key, "2"), track.getAttributeColor(key, "3"));
    }

    @Test
    public void testSessionRoundTrip() {
        track.setColorScale(VariantTrack.formatKey("GQ"), GQ_SCALE);
        track.setColorGenotypesBy("GQ");

        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertEquals("GQ", json.getString("genotypeColorBy"));

        VariantTrack restored = new VariantTrack();
        restored.unmarshalJSON(json);
        assertEquals("GQ", restored.getColorGenotypesBy());
        assertEquals(VariantTrack.ColorMode.FORMAT_ATTRIBUTE, restored.getGenotypeColorMode());
        Genotype genotype = ((Variant) variants.get(0)).getGenotype("S3");
        assertEquals(GQ_SCALE.getColor(99), restored.getGenotypeAttributeColor(genotype));

        // Colors chosen for a categorical field
        track.setColorGenotypesBy("CL");
        Color high = colorOf(0, "S2");
        json = new JSONObject();
        track.marshalJSON(json);
        restored = new VariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(high, restored.getGenotypeAttributeColor(((Variant) variants.get(0)).getGenotype("S2")));
    }

    /**
     * The mode without a field to color by has nothing to draw, so a session holding one reverts to genotype.
     */
    @Test
    public void testModeWithoutFieldRevertsToGenotype() {
        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        json.put("genotypeColorMode", "FORMAT_ATTRIBUTE");

        VariantTrack restored = new VariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(VariantTrack.ColorMode.GENOTYPE, restored.getGenotypeColorMode());
    }

    /**
     * The renderer fills the genotype cell with the field's color.
     */
    @Test
    public void testRenderGenotype() {
        track.setColorScale(VariantTrack.formatKey("GQ"), GQ_SCALE);
        track.setColorGenotypesBy("GQ");

        BufferedImage image = new BufferedImage(20, 20, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        Rectangle rect = new Rectangle(0, 0, 20, 20);
        RenderContext context = new RenderContext(null, graphics, null, rect, rect, rect);

        new VariantRenderer(track).renderGenotypeBandSNP((Variant) variants.get(0), context, rect, 0, 20, "S3",
                track.getGenotypeColorMode(), false);

        assertEquals(GQ_SCALE.getColor(99), new Color(image.getRGB(2, 10)));
        assertEquals(GQ_SCALE.getColor(99), new Color(image.getRGB(18, 10)));
    }

    // A scale is a table of colors, so its ends are a step off the colors it was given
    private static void assertClose(Color expected, Color actual) {
        assertEquals(expected.getRed(), actual.getRed(), 5);
        assertEquals(expected.getGreen(), actual.getGreen(), 5);
        assertEquals(expected.getBlue(), actual.getBlue(), 5);
    }

    private Color colorOf(int index, String sample) {
        return track.getGenotypeAttributeColor(((Variant) variants.get(index)).getGenotype(sample));
    }
}
