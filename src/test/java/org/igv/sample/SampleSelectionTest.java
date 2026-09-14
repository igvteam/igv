package org.igv.sample;

import org.igv.AbstractHeadlessTest;
import org.igv.track.AttributeManager;
import org.igv.track.TrackLoader;
import org.igv.ui.SampleSelectionDialog;
import org.igv.util.FilterElement;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.igv.variant.VariantTrack;
import org.json.JSONObject;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;

/**
 * Tests for filtering a track's samples by a list of IDs (issue #1215).
 */
public class SampleSelectionTest extends AbstractHeadlessTest {

    private static final String FILE_PATH = TestUtils.DATA_DIR + "vcf/multi_allele_freqs.vcf";

    private VariantTrack track;

    @Before
    public void loadTrack() throws Exception {
        TestUtils.createIndex(FILE_PATH);
        track = loadVariantTrack();
        assertEquals(15, track.sampleCount());
    }

    @After
    public void clearAttributes() {
        AttributeManager.getInstance().clearAllAttributes();
    }

    @Test
    public void testParseSampleIds() {
        assertEquals(Arrays.asList("CC-124", "D66", "sample one", "2137"),
                SampleSelectionDialog.parseSampleIds(" CC-124\r\nD66,sample one\t2137\n\nD66\n"));
        assertNull(SampleSelectionDialog.parseSampleIds(" \n,\t"));
    }

    /**
     * Selected samples are shown in the order of the list, as in igv.js.
     */
    @Test
    public void testSelection() {
        track.setSelectedSamples(Arrays.asList("D66", "CC-124", "2137"));
        assertEquals(Arrays.asList("D66", "CC-124", "2137"), visibleSamples());

        track.setSelectedSamples(null);
        assertEquals(15, track.sampleCount());
    }

    @Test
    public void testSortWithSelection() {
        track.setSelectedSamples(Arrays.asList("D66", "CC-124", "2137"));
        track.sortSamples(Comparator.naturalOrder());
        assertEquals(Arrays.asList("2137", "CC-124", "D66"), visibleSamples());
        assertEquals(Arrays.asList("2137", "CC-124", "D66"), track.getSelectedSamples());
    }

    /**
     * The dialog opens listing every sample, so accepting it unchanged must not save a selection.
     */
    @Test
    public void testSelectingAllSamplesClearsSelection() {
        track.setSelectedSamples(Arrays.asList("D66", "CC-124"));
        track.setSelectedSamples(track.getSampleNames());
        assertNull(track.getSelectedSamples());
        assertEquals(15, track.sampleCount());

        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertFalse(json.has("samples"));
    }

    /**
     * "samples" is only the ID filter -- a sort is saved as "sort", not as the sample list.
     */
    @Test
    public void testSortDoesNotWriteSamples() {
        track.sortSamplesByName(false);
        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertFalse(json.has("samples"));
    }

    /**
     * Samples must pass both the ID filter and the attribute filter to be shown.
     */
    @Test
    public void testIdAndAttributeFiltersCombine() {
        AttributeManager attributeManager = AttributeManager.getInstance();
        attributeManager.addAttribute("CC-124", "strain", "lab");
        attributeManager.addAttribute("CC-125", "strain", "lab");
        attributeManager.addAttribute("D66", "strain", "wild");
        SampleFilter labFilter = new SampleFilter(true,
                List.of(new FilterElement("strain", FilterElement.Operator.EQUAL, "lab")));

        track.setSelectedSamples(Arrays.asList("D66", "CC-124"));
        track.setSampleFilter(labFilter);
        assertEquals(List.of("CC-124"), visibleSamples());

        // Each filter stays on when the other is cleared
        track.setSelectedSamples(null);
        assertEquals(Arrays.asList("CC-124", "CC-125"), visibleSamples());

        track.setSelectedSamples(Arrays.asList("D66", "CC-124"));
        track.setSampleFilter(null);
        assertEquals(Arrays.asList("D66", "CC-124"), visibleSamples());
    }

    /**
     * The attribute filter's match mode survives a session round trip.
     */
    @Test
    public void testSampleFilterMatchModeRoundTrip() {
        for (boolean matchAll : new boolean[]{true, false}) {
            SampleFilter filter = new SampleFilter(matchAll,
                    List.of(new FilterElement("strain", FilterElement.Operator.EQUAL, "lab")));
            JSONObject json = filter.toJson();
            assertEquals(matchAll ? "all" : "any", json.getString("match"));
            assertEquals(matchAll, SampleFilter.fromJson(json).isMatchAll());
        }
    }

    @Test
    public void testSessionRoundTrip() {
        track.setSelectedSamples(Arrays.asList("cw15", "CC-125"));

        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertEquals(Arrays.asList("cw15", "CC-125"), json.getJSONArray("samples").toList());

        VariantTrack restored = loadVariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(Arrays.asList("cw15", "CC-125"), restored.getSelectedSamples());
        assertEquals(Arrays.asList("cw15", "CC-125"), visibleSamples(restored));
    }

    /**
     * Both filters are saved and restored together.
     */
    @Test
    public void testCombinedFiltersRoundTrip() {
        AttributeManager.getInstance().addAttribute("CC-124", "strain", "lab");
        track.setSelectedSamples(Arrays.asList("D66", "CC-124"));
        track.setSampleFilter(new SampleFilter(true,
                List.of(new FilterElement("strain", FilterElement.Operator.EQUAL, "lab"))));

        JSONObject json = new JSONObject();
        track.marshalJSON(json);

        VariantTrack restored = loadVariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(Arrays.asList("D66", "CC-124"), restored.getSelectedSamples());
        assertEquals(List.of("CC-124"), visibleSamples(restored));
    }

    /**
     * An igv.js track configuration's "samples" filters the track, and the rest of the samples are still
     * available to show again.  IDs not in the file are ignored.
     */
    @Test
    public void testIgvJsSamplesProperty() {
        JSONObject json = new JSONObject("{\"samples\": [\"S1D2\", \"not-a-sample\", \"2137\"]}");
        track.unmarshalJSON(json);
        assertEquals(Arrays.asList("S1D2", "2137"), visibleSamples());
        assertEquals(15, track.getSampleNames().size());

        track.setSelectedSamples(null);
        assertEquals(15, track.sampleCount());
    }

    private VariantTrack loadVariantTrack() {
        return (VariantTrack) (new TrackLoader()).load(new ResourceLocator(FILE_PATH), genome).get(0);
    }

    private List<String> visibleSamples() {
        return visibleSamples(track);
    }

    private static List<String> visibleSamples(VariantTrack t) {
        return t.getSampleGroups().stream().flatMap(g -> g.samples().stream()).collect(Collectors.toList());
    }
}
