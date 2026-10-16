package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.igv.AbstractHeadlessTest;
import org.igv.DirectoryManager;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ContinuousColorScale;
import org.igv.track.TrackLoader;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
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
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

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
                "CADD_PHRED\tContinuousColorScale;0.0;40.0;255,255,200;255,0,0",
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

    @Test
    public void testBuiltinScheme() {
        assertEquals(new Color(202, 0, 32), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(255, 33, 1), VariantColorSchemes.getColor("SVTYPE", "DEL"));
        assertEquals(new Color(55, 126, 184), VariantColorSchemes.getColor("VT", "SNP"));
        assertNull(VariantColorSchemes.getColor("CLNSIG", "drug_response"));
        assertNull(VariantColorSchemes.getColor("NOT_AN_ATTRIBUTE", "x"));
    }

    /**
     * A user's scheme for an attribute replaces the built-in one entirely -- there is one scheme per attribute.
     */
    @Test
    public void testUserSchemeReplacesBuiltin() throws Exception {
        writeScheme("CLNSIG.txt", "#name=Mine", "CLNSIG\tPathogenic\t1,2,3");
        VariantColorSchemes.reset();

        assertEquals(new Color(1, 2, 3), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        // Values the user scheme does not cover do not fall through to the built-in
        assertNull(VariantColorSchemes.getColor("CLNSIG", "Benign"));
        assertEquals("One scheme per attribute", 1,
                VariantColorSchemes.getSchemes().stream().filter(s -> "CLNSIG".equals(s.getKey())).count());
        assertFalse(VariantColorSchemes.getScheme("CLNSIG").isBuiltIn());
    }


    @Test
    public void testRemove() throws Exception {
        writeScheme("CLNSIG.txt", "#name=Mine", "CLNSIG\tPathogenic\t1,2,3");
        VariantColorSchemes.reset();

        VariantColorScheme scheme = VariantColorSchemes.getUserSchemes().get(0);
        assertTrue(VariantColorSchemes.remove(scheme));
        assertFalse(scheme.getFile().exists());
        assertTrue(VariantColorSchemes.getUserSchemes().isEmpty());

        // Back to the built-in color
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

        writeScheme("AF.txt", "AF\tContinuousColorScale;0.0;0.1;255,255,200;255,0,0");
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

        track.setAttributeColor("CLNSIG", "Pathogenic", new Color(7, 8, 9));
        assertEquals(new Color(7, 8, 9), track.getAttributeColor((Variant) variants.get(0)));

        track.clearAttributeColors("CLNSIG");
        assertEquals(new Color(202, 0, 32), track.getAttributeColor((Variant) variants.get(0)));
    }

    /**
     * Saving the legend writes a scheme that applies to every track, not just the one it was edited on.
     */
    @Test
    public void testSaveScheme() throws Exception {
        Map<String, Color> colors = new LinkedHashMap<>();
        colors.put("Pathogenic", new Color(1, 2, 3));
        colors.put("drug_response", new Color(4, 5, 6));

        VariantColorScheme scheme = VariantColorSchemes.saveScheme("CLNSIG", colors);

        assertEquals("CLNSIG.txt", scheme.getFile().getName());
        assertTrue(scheme.getFile().exists());

        // Shadows the built-in, and covers a value the built-in does not
        VariantColorSchemes.reset();
        assertEquals(new Color(1, 2, 3), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
        assertEquals(new Color(4, 5, 6), VariantColorSchemes.getColor("CLNSIG", "drug_response"));
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
     * A double gradient round trips exactly, including the neutral band between negStart and posStart, which
     * a three number range could not represent.
     */
    @Test
    public void testDoubleGradientScale() throws Exception {
        // (negStart, negEnd, posStart, posEnd) -- negEnd is the minimum, negStart the inner edge of the band
        ContinuousColorScale scale = new ContinuousColorScale(-2, -10, 2, 10,
                new Color(0, 0, 255), new Color(255, 255, 255), new Color(255, 0, 0));

        VariantColorSchemes.saveScale("SCORE", scale);
        VariantColorSchemes.reset();

        ContinuousColorScale restored = (ContinuousColorScale) VariantColorSchemes.getScale("SCORE");
        assertNotNull(restored);
        assertTrue(restored.isUseDoubleGradient());
        assertEquals(-10.0, restored.getMinimum(), 1e-9);
        assertEquals(-2.0, restored.getNegStart(), 1e-9);
        assertEquals(2.0, restored.getPosStart(), 1e-9);
        assertEquals(10.0, restored.getMaximum(), 1e-9);
        assertEquals(new Color(255, 255, 255), restored.getColor(0f));
    }

    /**
     * The answer to "is this numeric attribute categorical?" is stored in the scheme, so it is asked once.  It
     * has to be stored as a declaration, not merely as colors -- there may be no values loaded to color.
     */
    @Test
    public void testCategoricalDeclaration() throws Exception {
        assertFalse(VariantColorSchemes.getKeys().contains("DP"));

        VariantColorSchemes.saveScheme("DP", Collections.emptyMap(), true);
        VariantColorSchemes.reset();

        assertTrue(VariantColorSchemes.getUserSchemes().get(0).isCategorical("DP"));
        // The key counts as covered, which is what stops the scale dialog reappearing
        assertTrue(VariantColorSchemes.getKeys().contains("DP"));
        assertNull(VariantColorSchemes.getScale("DP"));
        // No colors were recorded, so values still take them from the palette
        assertNull(VariantColorSchemes.getColor("DP", "17"));

        // The declaration reads as plain text, so a scheme file can be written by hand
        VariantColorScheme scheme = VariantColorScheme.parse(
                new BufferedReader(new StringReader("TIER\tcategorical\nTIER\t1\t200,0,0\n")), "test");
        assertTrue(scheme.isCategorical("TIER"));
        assertEquals(new Color(200, 0, 0), scheme.getColor("TIER", "1"));
    }

    /**
     * "Number=A" and "Number=R" attributes carry one value per allele, so a multi-allelic record arrives as a
     * comma separated list.  Parsing the whole string as one number fails, which would render every such record
     * missing-gray even with a scale defined.  An allele frequency sums -- the total across alternate alleles is
     * what the attribute means, and it matches the "Allele Frequency" color mode.
     */
    @Test
    public void testAlleleFrequencyIsSummed() throws Exception {
        writeScheme("AF.txt", "AF\tContinuousColorScale;0.0;1.0;255,255,204;202,0,32");
        VariantColorSchemes.reset();

        VariantTrack track = loadTrack();
        track.setColorByAttribute("AF");

        Color expected = track.getAttributeColor("AF", "0.3");
        assertEquals(expected, track.getAttributeColor("AF", "0.1,0.2"));
        assertEquals(expected, track.getAttributeColor("AF", "0.2,0.1"));

        // A missing element does not make the whole record unusable
        assertEquals(track.getAttributeColor("AF", "0.2"), track.getAttributeColor("AF", "0.2,."));
        assertEquals(track.getAttributeColor("AF", "0.2"), track.getAttributeColor("AF", ".,0.2"));

        assertNotEquals(Color.gray, track.getAttributeColor("AF", "0.1,0.2"));
    }

    /**
     * Summing is only right for a frequency.  Anything else -- a score, a depth -- takes its maximum, so a
     * multi-allelic record never runs off the end of the scale.
     */
    @Test
    public void testOtherAttributesTakeTheMaximum() throws Exception {
        writeScheme("ALLELEID.txt", "ALLELEID\tContinuousColorScale;0.0;100.0;255,255,204;202,0,32");
        VariantColorSchemes.reset();

        VariantTrack track = loadTrack();
        track.setColorByAttribute("ALLELEID");

        assertEquals(track.getAttributeColor("ALLELEID", "40"), track.getAttributeColor("ALLELEID", "10,40"));
        assertEquals(track.getAttributeColor("ALLELEID", "40"), track.getAttributeColor("ALLELEID", "40,10"));
    }

    /**
     * Number=R lists the reference allele's value first.  It is not an alternate, so it must not be counted --
     * neither into a sum that claims to be the non-reference total, nor as the maximum.
     */
    @Test
    public void testReferenceValueIsSkipped() throws Exception {
        writeScheme("RDP.txt", "RDP\tContinuousColorScale;0.0;100.0;255,255,204;202,0,32");
        VariantColorSchemes.reset();

        VariantTrack track = loadTrack();
        track.setColorByAttribute("RDP");

        // 90 is the reference; the alternates are 5 and 20.  "0,20" is a record whose only alternate is 20, so it
        // colors the same -- and "0,90" does not, which shows the reference really is left out.
        assertEquals(track.getAttributeColor("RDP", "0,20"), track.getAttributeColor("RDP", "90,5,20"));
        assertNotEquals(track.getAttributeColor("RDP", "0,90"), track.getAttributeColor("RDP", "90,5,20"));

        // A record with no alternate allele (ALT is ".") has only the reference value -- nothing to color by.
        // Coloring by that value would make a reference-only site look like a variant with that score.
        assertEquals(Color.gray, track.getAttributeColor("RDP", "90"));
    }

    /**
     * A sum past the top of the scale saturates rather than falling back to the missing color.
     */
    @Test
    public void testSumBeyondScaleSaturates() throws Exception {
        writeScheme("AF.txt", "AF\tContinuousColorScale;0.0;1.0;255,255,204;202,0,32");
        VariantColorSchemes.reset();

        VariantTrack track = new VariantTrack();
        track.setColorByAttribute("AF");

        assertEquals(track.getAttributeColor("AF", "1.0"), track.getAttributeColor("AF", "0.7,0.8"));
    }

    /**
     * The IGV directory can move (Preferences > Advanced).  Schemes cached from the old one hold paths into a
     * directory that no longer exists, so saving would fail and removals would reappear.
     */
    @Test
    public void testCacheFollowsTheIgvDirectory() throws Exception {
        writeScheme("CLNSIG.txt", "#name=Mine", "CLNSIG\tPathogenic\t1,2,3");
        VariantColorSchemes.reset();
        assertEquals(new Color(1, 2, 3), VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));

        File moved = new File(TestUtils.TMP_OUTPUT_DIR, "igv-moved");
        moved.mkdirs();
        DirectoryManager.setIgvDirectory(moved);

        // No explicit reset -- the cache has to notice by itself
        assertTrue(VariantColorSchemes.getUserSchemes().isEmpty());
        assertEquals("Back to the built-in color", new Color(202, 0, 32),
                VariantColorSchemes.getColor("CLNSIG", "Pathogenic"));
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
