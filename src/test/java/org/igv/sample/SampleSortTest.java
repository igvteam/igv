package org.igv.sample;

import org.igv.AbstractHeadlessTest;
import org.igv.seg.SegTrack;
import org.igv.track.AbstractTrack;
import org.igv.track.AttributeManager;
import org.igv.track.RegionScoreType;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.igv.variant.Variant;
import org.igv.variant.VariantTrack;
import org.json.JSONObject;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Tests for saving and restoring sample sorts in the igv.js "sort" format.
 */
public class SampleSortTest extends AbstractHeadlessTest {

    private static final String VCF_PATH = TestUtils.DATA_DIR + "vcf/multi_allele_freqs.vcf";
    private static final String SEG_PATH = TestUtils.DATA_DIR + "seg/canFam2_hg18.seg";

    private VariantTrack track;
    private Variant variant;

    @Before
    public void loadTrack() throws Exception {
        TestUtils.createIndex(VCF_PATH);
        track = loadVariantTrack();
        variant = (Variant) track.getFeatures("chr1", 0, 1000000).get(2);    // chr1:543230, multi-allelic
    }

    @After
    public void clearAttributes() {
        AttributeManager.getInstance().clearAllAttributes();
    }

    @Test
    public void testIgvJsAttributeFormat() {
        SampleSort sort = SampleSort.fromJson(new JSONObject("{\"option\": \"ATTRIBUTE\", \"attribute\": \"Population\", \"direction\": \"ASC\"}"));
        assertArrayEquals(new String[]{"Population"}, sort.getAttributes());
        assertTrue(sort.getAttributeAscending()[0]);

        JSONObject json = sort.toJson();
        assertEquals("Population", json.getString("attribute"));
        assertEquals("ASC", json.getString("direction"));
    }

    /**
     * More than one attribute is written as parallel arrays.
     */
    @Test
    public void testMultiAttributeFormat() {
        SampleSort sort = SampleSort.attributes(new String[]{"Population", "Sex"}, new boolean[]{true, false});
        JSONObject json = sort.toJson();
        assertEquals(List.of("Population", "Sex"), json.getJSONArray("attribute").toList());
        assertEquals(List.of("ASC", "DESC"), json.getJSONArray("direction").toList());

        SampleSort restored = SampleSort.fromJson(json);
        assertArrayEquals(new String[]{"Population", "Sex"}, restored.getAttributes());
        assertTrue(restored.getAttributeAscending()[0]);
        assertFalse(restored.getAttributeAscending()[1]);
    }

    /**
     * igv.js writes the genotype option in lower case, and accepts a 1-based position in place of start and end.
     */
    @Test
    public void testIgvJsLocusFormats() {
        SampleSort sort = SampleSort.fromJson(new JSONObject("{\"option\": \"genotype\", \"direction\": \"DESC\", \"chr\": \"chr1\", \"start\": 10, \"end\": 20}"));
        assertEquals(SampleSort.GENOTYPE, sort.getOption());
        assertFalse(sort.isAscending());
        assertEquals(10, sort.getStart());
        assertEquals(20, sort.getEnd());

        sort = SampleSort.fromJson(new JSONObject("{\"direction\": \"ASC\", \"chr\": \"chr1\", \"position\": 100}"));
        assertNull(sort.getOption());
        assertTrue(sort.isAscending());
        assertEquals(99, sort.getStart());
        assertEquals(100, sort.getEnd());
    }

    /**
     * Only the last sort is saved, and a restored sort starts from file order, so each case starts from a fresh track --
     * samples tied in one sort would otherwise keep the order of the previous sort.
     */
    @Test
    public void testVariantSortsRoundTrip() {
        for (String option : new String[]{SampleSort.GENOTYPE, SampleSort.DEPTH, SampleSort.QUALITY}) {
            for (boolean ascending : new boolean[]{true, false}) {
                track = loadVariantTrack();
                variant = (Variant) track.getFeatures("chr1", 0, 1000000).get(2);
                track.sortSamples(option, variant, ascending);
                JSONObject json = new JSONObject();
                track.marshalJSON(json);

                JSONObject sortJson = json.getJSONObject("sort");
                assertEquals(option, sortJson.getString("option"));
                assertEquals(ascending ? "ASC" : "DESC", sortJson.getString("direction"));
                assertEquals(variant.getStart(), sortJson.getInt("start"));
                assertFalse(json.has("samples"));

                assertRestoredOrder(json, option + " " + sortJson.getString("direction"));
            }
        }

        // A sort by sample name, which needs no variant
        for (boolean ascending : new boolean[]{true, false}) {
            track = loadVariantTrack();
            track.sortSamplesByName(ascending);
            JSONObject json = new JSONObject();
            track.marshalJSON(json);
            assertEquals(SampleSort.SAMPLE_NAME, json.getJSONObject("sort").getString("option"));
            assertRestoredOrder(json, "SAMPLE_NAME " + ascending);
        }

        // A sort by several attributes, written as parallel arrays
        track = loadVariantTrack();
        AttributeManager attributeManager = AttributeManager.getInstance();
        for (String sample : track.getSampleNames()) {
            attributeManager.addAttribute(sample, "group", sample.startsWith("CC") ? "b" : "a");
            attributeManager.addAttribute(sample, "rank", String.valueOf(sample.length()));
        }
        track.sortSamplesByAttributes(new String[]{"group", "rank"}, new boolean[]{false, true});
        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertEquals(List.of("group", "rank"), json.getJSONObject("sort").getJSONArray("attribute").toList());
        assertTrue(visibleSamples(track).get(0).startsWith("CC"));
        assertRestoredOrder(json, "ATTRIBUTE");
    }

    /**
     * The genotype sort at the test variant changes the file order, so the round-trip tests are not trivially true.
     */
    @Test
    public void testGenotypeSortChangesOrder() {
        List<String> fileOrder = visibleSamples(track);
        track.sortSamples(SampleSort.GENOTYPE, variant, false);
        assertNotEquals(fileOrder, visibleSamples(track));
    }

    /**
     * A genotype sort written by igv.js restores the same order as the desktop sort in that direction.
     */
    @Test
    public void testIgvJsGenotypeSort() {
        track.sortSamples(SampleSort.GENOTYPE, variant, false);
        JSONObject json = new JSONObject("{\"sort\": {\"option\": \"genotype\", \"direction\": \"DESC\", \"chr\": \"chr1\", " +
                "\"start\": " + (variant.getStart() - 2) + ", \"end\": " + (variant.getEnd() + 2) + "}}");
        VariantTrack restored = loadVariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(visibleSamples(track), visibleSamples(restored));
    }

    /**
     * Seg region sorts are saved as igv.js VALUE sorts -- deletion is ascending, the others descending -- and
     * survive a later ID filter, so the session restores the order shown.
     */
    @Test
    public void testSegValueSortRoundTrip() {
        SegTrack segTrack = loadSegTrack();
        for (RegionScoreType type : new RegionScoreType[]{RegionScoreType.AMPLIFICATION, RegionScoreType.DELETION}) {
            segTrack.sortSamplesByValue("chr1", 0, 250000000, type);
            JSONObject json = new JSONObject();
            segTrack.marshalJSON(json);

            JSONObject sortJson = json.getJSONObject("sort");
            assertEquals(SampleSort.VALUE, sortJson.getString("option"));
            assertEquals(type == RegionScoreType.DELETION ? "ASC" : "DESC", sortJson.getString("direction"));

            SegTrack restored = loadSegTrack();
            restored.unmarshalJSON(json);
            assertEquals(type.toString(), visibleSamples(segTrack), visibleSamples(restored));
        }

        // An ID filter applied after the sort keeps the sorted order, and that is what the session restores
        List<String> sorted = visibleSamples(segTrack);
        List<String> subset = Arrays.asList(sorted.get(3), sorted.get(0), sorted.get(2));
        segTrack.setSelectedSamples(subset);
        List<String> expected = sorted.stream().filter(subset::contains).collect(Collectors.toList());
        assertEquals(expected, visibleSamples(segTrack));

        JSONObject json = new JSONObject();
        segTrack.marshalJSON(json);
        SegTrack restored = loadSegTrack();
        restored.unmarshalJSON(json);
        assertEquals(expected, visibleSamples(restored));
    }

    /**
     * An ID filter and a sort combine whichever order they are applied in: the filtered samples are shown in sort
     * order, and the session restores the order shown.
     */
    @Test
    public void testVariantSortAndIdFilter() {
        track.sortSamples(SampleSort.GENOTYPE, variant, false);
        List<String> sorted = visibleSamples(track);

        List<String> subset = Arrays.asList(sorted.get(5), sorted.get(1), sorted.get(10), sorted.get(0));
        track.setSelectedSamples(subset);
        assertEquals(sorted.stream().filter(subset::contains).collect(Collectors.toList()), visibleSamples(track));
        assertRestoredOrder(sessionJson(), "sort then ID filter");

        // The other order -- the filter first, then a sort over the filtered samples
        track = loadVariantTrack();
        track.setSelectedSamples(Arrays.asList("D66", "CC-124", "2137", "cw15"));
        track.sortSamplesByName(true);
        JSONObject json = sessionJson();
        assertEquals(List.of("2137", "CC-124", "D66", "cw15"), json.getJSONArray("samples").toList());
        assertRestoredOrder(json, "ID filter then sort");
    }

    /**
     * A sort that can't be read is skipped; the rest of the track's settings are restored.
     */
    @Test
    public void testUnreadableSortIsSkipped() {
        JSONObject json = new JSONObject("{\"samples\": [\"D66\", \"2137\"], \"sort\": {\"option\": \"GENOTYPE\", \"direction\": \"ASC\"}}");
        track.unmarshalJSON(json);
        assertNull(track.getSampleSort());
        assertEquals(List.of("D66", "2137"), visibleSamples(track));
    }

    private JSONObject sessionJson() {
        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        return json;
    }

    private void assertRestoredOrder(JSONObject json, String message) {
        VariantTrack restored = loadVariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(message, visibleSamples(track), visibleSamples(restored));
    }

    private VariantTrack loadVariantTrack() {
        return (VariantTrack) (new TrackLoader()).load(new ResourceLocator(VCF_PATH), genome).get(0);
    }

    private SegTrack loadSegTrack() {
        return (SegTrack) (new TrackLoader()).load(new ResourceLocator(SEG_PATH), genome).stream()
                .filter(t -> t instanceof SegTrack).findFirst().orElseThrow();
    }

    private static List<String> visibleSamples(AbstractTrack t) {
        return t.getSampleGroups().stream().flatMap(g -> g.samples().stream()).collect(Collectors.toList());
    }
}
