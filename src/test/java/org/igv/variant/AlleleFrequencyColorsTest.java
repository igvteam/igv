package org.igv.variant;

import htsjdk.tribble.Feature;
import org.igv.AbstractHeadlessTest;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.ColorStopScale;
import org.igv.track.RenderContext;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.Test;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests for coloring variants by allele frequency rarity.
 */
public class AlleleFrequencyColorsTest extends AbstractHeadlessTest {

    private static final String SITES = TestUtils.DATA_DIR + "vcf/allele_frequency_sites.vcf";
    private static final String GENOTYPES = TestUtils.DATA_DIR + "vcf/allele_frequency_genotypes.vcf";

    private static final ColorStopScale SCALE = AlleleFrequencyColors.DEFAULT_SCALE;
    private static final Color MOST_RARE = SCALE.getStops().get(0).color();

    /**
     * AF, then GMAF, then annotation fields such as gnomAD_AF; the rarest allele at a multi-allelic site.
     */
    @Test
    public void testFrequencies() {
        List<Variant> variants = variants(load(SITES));
        double[] expected = {0.2, 0.0002, 0.004, 0.00005, -1};
        for (int i = 0; i < expected.length; i++) {
            assertEquals("variant " + i, expected[i],
                    AlleleFrequencyColors.getFrequency(variants.get(i), VariantTrack.ColorMode.ALLELE_FREQUENCY), 1e-12);
        }
    }

    @Test
    public void testAlleleFractions() {
        List<Variant> variants = variants(load(GENOTYPES));
        double[] expected = {0.25, 0.0005, -1};
        for (int i = 0; i < expected.length; i++) {
            assertEquals("variant " + i, expected[i],
                    AlleleFrequencyColors.getFrequency(variants.get(i), VariantTrack.ColorMode.ALLELE_FRACTION), 1e-12);
        }
    }

    /**
     * Common variants are neutral, rarer ones shift color, and a variant with no frequency is colored as the rarest.
     */
    @Test
    public void testColors() {
        VariantTrack track = load(SITES);
        track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
        List<Variant> variants = variants(track);

        Color common = track.getAlleleFrequencyColor(variants.get(0));
        assertEquals(SCALE.getStops().get(SCALE.getStops().size() - 1).color(), common);
        assertEquals(SCALE.getColor(0.0002), track.getAlleleFrequencyColor(variants.get(1)));
        assertNotEquals(common, track.getAlleleFrequencyColor(variants.get(1)));
        assertEquals(MOST_RARE, track.getAlleleFrequencyColor(variants.get(4)));
    }

    /**
     * A file that declares no field a mode reads isn't colored as rare -- its variants take the track color.
     */
    @Test
    public void testNoFrequencyFields() {
        VariantTrack track = load(GENOTYPES);
        List<Variant> variants = variants(track);

        track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
        assertFalse(track.hasFrequencyFields(VariantTrack.ColorMode.ALLELE_FREQUENCY));
        assertEquals(track.getColor(), track.getAlleleFrequencyColor(variants.get(0)));

        track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FRACTION);
        assertTrue(track.hasFrequencyFields(VariantTrack.ColorMode.ALLELE_FRACTION));
        assertEquals(SCALE.getColor(0.25), track.getAlleleFrequencyColor(variants.get(0)));
        assertEquals(MOST_RARE, track.getAlleleFrequencyColor(variants.get(2)));
    }

    /**
     * Rarity colors are the default, except for a file without genotypes that has allele frequencies.
     */
    @Test
    public void testDefaultDisplay() {
        assertEquals(AlleleFrequencyColors.DISPLAY_AUTOMATIC,
                PreferencesManager.getPreferences().get(Constants.VARIANT_ALLELE_FREQUENCY_DISPLAY));
        assertTrue(load(SITES).isAlleleFrequencyBars());
        assertFalse(load(GENOTYPES).isAlleleFrequencyBars());
    }

    /**
     * Whether AF describes the file's own samples can't be known, so the preference can choose either display for
     * every file, and a choice made on a track overrides it.
     */
    @Test
    public void testDisplayPreference() {
        IGVPreferences prefs = PreferencesManager.getPreferences();
        try {
            prefs.put(Constants.VARIANT_ALLELE_FREQUENCY_DISPLAY, AlleleFrequencyColors.DISPLAY_COLOR_SCALE);
            assertFalse(load(SITES).isAlleleFrequencyBars());
            assertFalse(load(GENOTYPES).isAlleleFrequencyBars());

            prefs.put(Constants.VARIANT_ALLELE_FREQUENCY_DISPLAY, AlleleFrequencyColors.DISPLAY_BAR);
            assertTrue(load(SITES).isAlleleFrequencyBars());
            assertTrue(load(GENOTYPES).isAlleleFrequencyBars());

            VariantTrack chosen = load(GENOTYPES);
            chosen.setAlleleFrequencyBars(false);
            assertFalse(chosen.isAlleleFrequencyBars());
        } finally {
            prefs.remove(Constants.VARIANT_ALLELE_FREQUENCY_DISPLAY);
        }
    }

    /**
     * Colored by rarity, the whole variant band takes the rarity color.
     */
    @Test
    public void testRenderRarityColors() {
        VariantTrack track = load(SITES);
        track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
        track.setAlleleFrequencyBars(false);
        Variant rare = variants(track).get(1);

        BufferedImage image = renderSiteBand(track, rare);
        Color expected = track.getAlleleFrequencyColor(rare);
        assertEquals(expected, new Color(image.getRGB(10, 4)));
        assertEquals(expected, new Color(image.getRGB(10, 24)));
    }

    /**
     * The bar display is unchanged: the band is split between the variant color, below, and the reference color.
     */
    @Test
    public void testRenderBars() {
        VariantTrack track = load(SITES);
        track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
        track.setAlleleFrequencyBars(true);

        BufferedImage image = renderSiteBand(track, variants(track).get(0));    // AF = 0.2
        assertNotEquals(new Color(image.getRGB(10, 4)), new Color(image.getRGB(10, 24)));
    }

    private static BufferedImage renderSiteBand(VariantTrack track, Variant variant) {
        BufferedImage image = new BufferedImage(20, 25, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        Rectangle rect = new Rectangle(0, 0, 20, 25);
        RenderContext context = new RenderContext(null, graphics, null, rect, rect, rect);
        new VariantRenderer(track).renderSiteBand(variant, rect, 0, 20, context);
        return image;
    }

    @Test
    public void testSessionRoundTrip() {
        VariantTrack track = load(SITES);
        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertFalse("Not written unless chosen for the track", json.has("alleleFrequencyBars"));

        track.setAlleleFrequencyBars(false);
        json = new JSONObject();
        track.marshalJSON(json);
        assertFalse(json.getBoolean("alleleFrequencyBars"));

        VariantTrack restored = load(SITES);
        restored.unmarshalJSON(json);
        assertFalse(restored.isAlleleFrequencyBars());
    }

    private VariantTrack load(String path) {
        return (VariantTrack) new TrackLoader().load(new ResourceLocator(path), genome).get(0);
    }

    private static List<Variant> variants(VariantTrack track) {
        List<Feature> features = track.getFeatures("chr1", 0, 1000);
        return features.stream().map(f -> (Variant) f).toList();
    }
}
