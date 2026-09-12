package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.igv.AbstractHeadlessTest;
import org.igv.DirectoryManager;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ContinuousColorScale;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Tests for variant color schemes -- the files that assign colors to VCF INFO attribute values.
 */
public class VariantColorSchemeTest extends AbstractHeadlessTest {

    private File previousIgvDirectory;
    private File igvDirectory;

    @Before
    public void useTemporaryIgvDirectory() throws Exception {
        previousIgvDirectory = DirectoryManager.getIgvDirectory();
        igvDirectory = new File(TestUtils.TMP_OUTPUT_DIR, "igv");
        igvDirectory.mkdirs();
        DirectoryManager.setIgvDirectory(igvDirectory);
        VariantColorSchemes.reset();
    }

    @After
    public void restoreIgvDirectory() {
        DirectoryManager.setIgvDirectory(previousIgvDirectory);
        VariantColorSchemes.reset();
    }

    @Test
    public void testParse() throws Exception {
        String contents = String.join("\n",
                "#name=Test scheme",
                "#description=A scheme for testing",
                "# a comment",
                "#colors",
                "CLNSIG\tPathogenic\t255,0,0",
                "CLNSIG\t*\t10,10,10",
                "CADD_PHRED\t0:40\t255,255,200\t255,0,0",
                "");

        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader(contents)), "fallback");

        assertEquals("Test scheme", scheme.getName());
        assertEquals("A scheme for testing", scheme.getDescription());
        assertEquals(new Color(255, 0, 0), scheme.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(255, 0, 0), scheme.getColor("clnsig", "pathogenic"));   // case insensitive
        assertEquals(new Color(10, 10, 10), scheme.getColor("CLNSIG", "Benign"));      // wildcard default
        assertNull(scheme.getColor("SVTYPE", "DEL"));                                  // key not covered
        assertEquals(Set.of("CLNSIG", "CADD_PHRED"), Set.copyOf(scheme.getKeys()));
    }

    /**
     * One malformed row does not reject the whole scheme.
     */
    @Test
    public void testParseSkipsBadRows() throws Exception {
        String contents = String.join("\n",
                "CLNSIG\tPathogenic\tnot-a-color",
                "CLNSIG\ttoo-few-fields",
                "CLNSIG\tBenign\t0,0,255",
                "");

        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader(contents)), "test");

        assertNull(scheme.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(0, 0, 255), scheme.getColor("CLNSIG", "Benign"));
    }

    /**
     * The name falls back to the file name when the file has no "#name=" directive.
     */
    @Test
    public void testDefaultName() throws Exception {
        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader("")), "lab-tiers");
        assertEquals("lab-tiers", scheme.getName());
    }

    @Test
    public void testBuiltinScheme() {
        assertEquals(new Color(202, 0, 32), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(255, 33, 1), VariantColorSchemes.getColor("SVTYPE", "DEL"));
        assertEquals(new Color(55, 126, 184), VariantColorSchemes.getColor("VT", "SNP"));
        assertNull(VariantColorSchemes.getColor("CLNSIG", "drug_response"));
        assertNull(VariantColorSchemes.getColor("NOT_AN_ATTRIBUTE", "x"));
    }

    /**
     * An imported scheme is searched before the built-in one, so a user can override IGV's colors.
     */
    @Test
    public void testUserSchemeShadowsBuiltin() throws Exception {
        writeScheme("mine.txt", "#name=Mine", "CLNSIG\tPathogenic\t1,2,3");
        VariantColorSchemes.reset();

        assertEquals(new Color(1, 2, 3), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        // Values the user scheme does not cover still come from the built-in
        assertEquals(new Color(5, 113, 176), VariantColorSchemes.getColor("CLNSIG", "Benign"));
    }

    /**
     * Importing copies the file into the IGV directory, so the scheme survives the original being deleted and is
     * there on the next startup.
     */
    @Test
    public void testImportCopiesFile() throws Exception {
        File source = new File(TestUtils.TMP_OUTPUT_DIR, "shared-colors.txt");
        try (PrintWriter writer = new PrintWriter(source)) {
            writer.println("#name=Shared");
            writer.println("TIER\t1\t200,0,0");
        }

        VariantColorScheme scheme = VariantColorSchemes.importFile(source);
        assertEquals("Shared", scheme.getName());
        assertEquals(source.getAbsolutePath(), scheme.getSource());

        File copy = new File(VariantColorSchemes.getSchemeDirectory(), "shared-colors.txt");
        assertTrue("Scheme was not copied into the IGV directory", copy.exists());

        // Delete the original -- the scheme is still there after a restart
        Files.delete(source.toPath());
        VariantColorSchemes.reset();
        assertEquals(new Color(200, 0, 0), VariantColorSchemes.getColor("TIER", "1"));
        assertEquals(List.of("Shared"),
                VariantColorSchemes.getUserSchemes().stream().map(VariantColorScheme::getName).collect(Collectors.toList()));
    }

    /**
     * Reading schemes must not create the scheme directory -- opening preferences should not leave a directory
     * behind for a feature the user never used.
     */
    @Test
    public void testReadingDoesNotCreateDirectory() {
        assertFalse(VariantColorSchemes.getSchemeDirectory().exists());
        VariantColorSchemes.getColor("CLNSIG", "Pathogenic");
        VariantColorSchemes.getSchemes();
        assertFalse("Reading schemes created the scheme directory", VariantColorSchemes.getSchemeDirectory().exists());
    }

    @Test
    public void testRemove() throws Exception {
        writeScheme("mine.txt", "#name=Mine", "CLNSIG\tPathogenic\t1,2,3");
        VariantColorSchemes.reset();

        VariantColorScheme scheme = VariantColorSchemes.getUserSchemes().get(0);
        assertTrue(VariantColorSchemes.remove(scheme));
        assertFalse(scheme.getFile().exists());
        assertTrue(VariantColorSchemes.getUserSchemes().isEmpty());

        // Back to the built-in color
        assertEquals(new Color(202, 0, 32), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
    }

    /**
     * Schemes shipped with IGV cannot be removed.
     */
    @Test
    public void testBuiltinCannotBeRemoved() {
        VariantColorScheme builtin = VariantColorSchemes.getBuiltinSchemes().get(0);
        assertTrue(builtin.isBuiltIn());
        assertFalse(VariantColorSchemes.remove(builtin));
        assertEquals(new Color(202, 0, 32), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
    }

    /**
     * A numeric attribute shades across the range a scheme gives it.
     */
    @Test
    public void testContinuousScale() throws Exception {
        String filePath = TestUtils.DATA_DIR + "vcf/clinvar_info.vcf";
        TestUtils.createIndex(filePath);
        VariantTrack track = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(filePath), genome).get(0);
        List<Feature> variants = track.getFeatures("chr1", 0, 1000);

        // AF is offered, but has no scale until a scheme gives it one -- selecting it asks for one
        assertTrue(colorableIds(track).contains("AF"));
        assertNull(VariantColorSchemes.getScale("AF"));

        writeScheme("af.txt", "AF\t0:0.1\t255,255,200\t255,0,0");
        VariantColorSchemes.reset();

        assertNotNull(VariantColorSchemes.getScale("AF"));
        track.setColorByAttribute("AF");

        Color low = track.getAttributeColor((Variant) variants.get(0));   // AF=0.01
        Color high = track.getAttributeColor((Variant) variants.get(8));  // AF=0.09
        assertNotEquals(low, high);
        assertTrue("Expected the high end of the scale to be redder", high.getRed() - high.getBlue() > low.getRed() - low.getBlue());
    }

    /**
     * Colors chosen on a track win over a scheme, and over palette colors assigned earlier.
     */
    @Test
    public void testTrackOverrideWinsOverScheme() throws Exception {
        VariantTrack track = loadTrack();
        List<Feature> variants = track.getFeatures("chr1", 0, 1000);
        track.setColorByAttribute("CLNSIG");

        assertEquals(new Color(202, 0, 32), track.getAttributeColor((Variant) variants.get(0)));   // built in

        track.setAttributeColorOverride("CLNSIG", "Pathogenic", new Color(7, 8, 9));
        assertEquals(new Color(7, 8, 9), track.getAttributeColor((Variant) variants.get(0)));

        track.clearAttributeColorOverrides("CLNSIG");
        assertEquals(new Color(202, 0, 32), track.getAttributeColor((Variant) variants.get(0)));
    }

    /**
     * The legend offers the values present in the loaded features.
     */
    @Test
    public void testAttributeValues() throws Exception {
        VariantTrack track = loadTrack();
        track.getFeatures("chr1", 0, 1000);
        track.setColorByAttribute("CLNSIG");

        // getFeatures does not pack features into the render cache, so nothing is "in view" until it does
        assertTrue(track.getAttributeValues("CLNSIG").isEmpty());

        track.setAttributeColorOverride("CLNSIG", "Pathogenic", Color.red);
        assertTrue(track.getAttributeColorOverrides("CLNSIG").containsKey("pathogenic"));
    }

    /**
     * Saving the legend writes a scheme that applies to every track, not just the one it was edited on.
     */
    @Test
    public void testSaveScheme() throws Exception {
        Map<String, Color> colors = new LinkedHashMap<>();
        colors.put("Pathogenic", new Color(1, 2, 3));
        colors.put("drug_response", new Color(4, 5, 6));

        VariantColorScheme scheme = VariantColorSchemes.saveScheme("My colors", "CLNSIG", colors);

        assertEquals("My colors", scheme.getName());
        assertTrue(scheme.getFile().exists());

        // Shadows the built-in, and covers a value the built-in does not
        VariantColorSchemes.reset();
        assertEquals(new Color(1, 2, 3), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(4, 5, 6), VariantColorSchemes.getColor("CLNSIG", "drug_response"));
    }

    /**
     * Colors the user chose survive a session round trip, separately from colors IGV assigned.
     */
    @Test
    public void testOverridesRoundTrip() throws Exception {
        VariantTrack track = loadTrack();
        track.setColorByAttribute("CLNSIG");
        track.setAttributeColorOverride("CLNSIG", "Pathogenic", new Color(7, 8, 9));

        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertEquals("7,8,9", json.getJSONObject("colorTable").getString("pathogenic"));

        VariantTrack restored = new VariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(new Color(7, 8, 9), restored.getAttributeColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(5, 113, 176), restored.getAttributeColor("CLNSIG", "Benign/Likely_benign"));
    }

    /**
     * A numeric attribute is a quantity, so it needs a color scale rather than a color per value.
     */
    @Test
    public void testNumericAttributes() throws Exception {
        VariantTrack track = loadTrack();
        track.getFeatures("chr1", 0, 1000);

        assertTrue(track.isNumericAttribute("AF"));           // Float
        assertTrue(track.isNumericAttribute("ALLELEID"));     // Integer
        assertFalse(track.isNumericAttribute("CLNSIG"));      // String
        assertFalse(track.isNumericAttribute("NOT_AN_ATTRIBUTE"));
    }

    /**
     * A saved scale round trips through the scheme file, in a form that can be edited by hand.
     */
    @Test
    public void testSaveScale() throws Exception {
        ContinuousColorScale scale =
                new ContinuousColorScale(0, 40, new Color(255, 255, 204), new Color(202, 0, 32));

        VariantColorSchemes.saveScale("CADD scale", "CADD_PHRED", scale);
        VariantColorSchemes.reset();

        AbstractColorScale restored = VariantColorSchemes.getScale("CADD_PHRED");
        assertNotNull(restored);
        assertEquals(scale.getColor(0f), restored.getColor(0f));
        assertEquals(scale.getColor(40f), restored.getColor(40f));
        assertNotEquals(restored.getColor(0f), restored.getColor(40f));

        assertEquals("CADD scale", VariantColorSchemes.getSchemeForScale("CADD_PHRED").getName());
    }

    /**
     * A three stop gradient, written as "min:mid:max" with three colors.
     */
    @Test
    public void testDoubleGradientScale() throws Exception {
        String contents = "SCORE\t-10:0:10\t0,0,255\t255,255,255\t255,0,0\n";
        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader(contents)), "test");

        AbstractColorScale scale = scheme.getScale("SCORE");
        assertNotNull(scale);
        assertNotEquals(scale.getColor(-10f), scale.getColor(10f));
        assertEquals(new Color(255, 255, 255), scale.getColor(0f));
    }

    /**
     * A three stop range needs three colors -- two is ambiguous, so the row is skipped rather than guessed at.
     */
    @Test
    public void testDoubleGradientNeedsThreeColors() throws Exception {
        String contents = "SCORE\t-10:0:10\t0,0,255\t255,0,0\n";
        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader(contents)), "test");
        assertNull(scheme.getScale("SCORE"));
    }

    /**
     * The answer to "is this numeric attribute categorical?" is stored in the scheme, so it is asked once.  It
     * has to be stored as a declaration, not merely as colors -- there may be no values in view to color.
     */
    @Test
    public void testCategoricalDeclaration() throws Exception {
        assertFalse(VariantColorSchemes.isCategorical("DP"));

        VariantColorSchemes.saveScheme("DP colors", "DP", Collections.emptyMap(), true);
        VariantColorSchemes.reset();

        assertTrue(VariantColorSchemes.isCategorical("DP"));
        // The key counts as covered, which is what stops the scale dialog reappearing
        assertTrue(VariantColorSchemes.getKeys().contains("DP"));
        assertNull(VariantColorSchemes.getScale("DP"));
        // No colors were recorded, so values still take them from the palette
        assertNull(VariantColorSchemes.getColor("DP", "17"));
    }

    /**
     * The declaration reads as plain text, so a scheme file can be written by hand.
     */
    @Test
    public void testCategoricalDeclarationParsed() throws Exception {
        String contents = "TIER\tcategorical\nTIER\t1\t200,0,0\n";
        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader(contents)), "test");

        assertTrue(scheme.isCategorical("TIER"));
        assertTrue(scheme.getKeys().contains("TIER"));
        assertEquals(new Color(200, 0, 0), scheme.getColor("TIER", "1"));
    }

    /**
     * A whole scheme round trips through the file, every attribute it covers -- the editor saves it this way.
     */
    @Test
    public void testSchemeRoundTrip() throws Exception {
        String contents = String.join("\n",
                "#name=Mixed",
                "#description=Several attributes at once",
                "#colors",
                "CLNSIG\tPathogenic\t1,2,3",
                "CLNSIG\t*\t9,9,9",
                "DP\tcategorical",
                "CADD\t0:40\t255,255,204\t202,0,32",
                "");
        VariantColorScheme scheme = VariantColorScheme.parse(new BufferedReader(new StringReader(contents)), "Mixed");

        VariantColorSchemes.save(scheme);
        VariantColorSchemes.reset();

        VariantColorScheme reloaded = VariantColorSchemes.getUserSchemes().get(0);
        assertEquals("Mixed", reloaded.getName());
        assertEquals("Several attributes at once", reloaded.getDescription());
        assertEquals(new Color(1, 2, 3), reloaded.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(9, 9, 9), reloaded.getColor("CLNSIG", "anything else"));
        assertTrue(reloaded.isCategorical("DP"));
        assertNotNull(reloaded.getScale("CADD"));
    }

    /**
     * A scheme shipped with IGV cannot be written over, so editing one saves a copy that shadows it.
     */
    @Test
    public void testEditingBuiltinSavesACopy() throws Exception {
        VariantColorScheme builtin = VariantColorSchemes.getBuiltinSchemes().get(0);
        assertNull(builtin.getFile());

        builtin.setColor("CLNSIG", "Pathogenic", new Color(1, 1, 1));
        VariantColorSchemes.save(builtin);
        VariantColorSchemes.reset();

        assertEquals(new Color(1, 1, 1), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        assertEquals(1, VariantColorSchemes.getUserSchemes().size());
        assertFalse(VariantColorSchemes.getUserSchemes().get(0).isBuiltIn());
    }

    private VariantTrack loadTrack() throws Exception {
        String filePath = TestUtils.DATA_DIR + "vcf/clinvar_info.vcf";
        TestUtils.createIndex(filePath);
        return (VariantTrack) (new TrackLoader()).load(new ResourceLocator(filePath), genome).get(0);
    }

    private List<String> colorableIds(VariantTrack track) {
        return track.getColorableInfoFields().stream().map(VCFInfoHeaderLine::getID).collect(Collectors.toList());
    }

    private void writeScheme(String fileName, String... lines) throws Exception {
        File directory = VariantColorSchemes.getSchemeDirectory();
        directory.mkdirs();
        File file = new File(directory, fileName);
        try (PrintWriter writer = new PrintWriter(file)) {
            for (String line : lines) {
                writer.println(line);
            }
        }
    }
}
