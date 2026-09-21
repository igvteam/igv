package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.igv.AbstractHeadlessTest;
import org.igv.DirectoryManager;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.ColorStopScale;
import org.igv.track.RenderContext;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Tests for coloring variants by allele frequency rarity.
 */
public class AlleleFrequencyColorsTest extends AbstractHeadlessTest {

    private static final String SITES = TestUtils.DATA_DIR + "vcf/allele_frequency_sites.vcf";
    private static final String GENOTYPES = TestUtils.DATA_DIR + "vcf/allele_frequency_genotypes.vcf";

    private static final ColorStopScale SCALE = AlleleFrequencyColors.DEFAULT_SCALE;
    private static final Color MOST_RARE = SCALE.getStops().get(0).color();

    private File previousIgvDirectory;

    /**
     * Color schemes the developer running the tests has installed must not cover the frequency fields tested here.
     */
    @Before
    public void isolateSchemes() {
        previousIgvDirectory = DirectoryManager.getIgvDirectory();
        File igvDirectory = new File(TestUtils.TMP_OUTPUT_DIR, "igv");
        igvDirectory.mkdirs();
        DirectoryManager.setIgvDirectory(igvDirectory);
        VariantColorSchemes.reset();
    }

    @After
    public void restoreIgvDirectory() {
        DirectoryManager.setIgvDirectory(previousIgvDirectory);
        VariantColorSchemes.reset();
    }

    /**
     * The top-level Allele Frequency reads AF, or GMAF if there is no AF, as the bar display does -- and no other
     * field.  The rarest allele at a multi-allelic site.
     */
    @Test
    public void testFrequencies() {
        List<Variant> variants = variants(load(SITES));
        // gnomAD_AF only, AF=., no fields, and the reference block have no value; the last has only GMAF
        double[] expected = {0.2, 0.0002, -1, -1, -1, -1, 0.03, -1};
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
     * AF has its own menu item, so it isn't offered again under color by INFO field.
     */
    @Test
    public void testMenuLeavesOutAF() {
        VariantTrack track = load(SITES);
        assertEquals(List.of("GMAF", "gnomAD_AF"),
                VariantTrackMenuHelper.getMenuInfoFields(track).stream().map(VCFInfoHeaderLine::getID).toList());
    }

    /**
     * Well-known allele frequency INFO fields are colored by rarity, without asking for a scale.
     */
    @Test
    public void testFrequencyInfoFields() {
        assertTrue(VariantColorSchemes.getScale("GMAF") instanceof ColorStopScale);
        assertTrue(VariantColorSchemes.getScale("gnomad_af") instanceof ColorStopScale);
        assertNull(VariantColorSchemes.getScale("AF"));

        VariantTrack track = load(SITES);
        assertTrue(VariantTrackMenuHelper.defineScaleIfNeeded(track, "gnomAD_AF"));

        track.setColorByAttribute("gnomAD_AF");
        List<Variant> variants = variants(track);
        assertEquals(AlleleFrequencyColors.getScale().getColor((float) 0.004),
                track.getAttributeColor(variants.get(2)));                          // gnomAD_AF=0.004

        // A multi-allelic site is colored by its rarest allele
        assertEquals(AlleleFrequencyColors.getScale().getColor((float) 0.00001),
                track.getAttributeColor(variants.get(7)));                          // gnomAD_AF=0.3,0.00001
    }

    /**
     * Rarity colors are the default, except for a file without genotypes that declares AF.
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
     * A reference block (only a <NON_REF> alternate) has no frequency but isn't a variant, so it keeps the non-ref
     * color rather than being colored as the rarest.
     */
    @Test
    public void testReferenceBlockColor() {
        VariantTrack track = load(SITES);
        track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
        track.setAlleleFrequencyBars(false);
        Variant referenceBlock = variants(track).get(5);
        assertTrue(referenceBlock.isNonRef());

        BufferedImage image = renderSiteBand(track, referenceBlock);
        assertEquals(new Color(200, 200, 215), new Color(image.getRGB(10, 10)));
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

    private static BufferedImage renderSiteBand(VariantTrack track, Variant variant) {
        BufferedImage image = new BufferedImage(20, 25, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        Rectangle rect = new Rectangle(0, 0, 20, 25);
        RenderContext context = new RenderContext(null, graphics, null, rect, rect, rect);
        new VariantRenderer(track).renderSiteBand(variant, rect, 0, 20, context);
        return image;
    }

    private VariantTrack load(String path) {
        return (VariantTrack) new TrackLoader().load(new ResourceLocator(path), genome).get(0);
    }

    private static List<Variant> variants(VariantTrack track) {
        List<Feature> features = track.getFeatures("chr1", 0, 1000);
        return features.stream().map(f -> (Variant) f).toList();
    }
}
