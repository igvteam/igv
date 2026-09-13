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
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

/**
 * Tests for restricting a track's samples to a list of IDs (issue #1215).
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
     * Selected samples are shown in track order, not the order they were entered.
     */
    @Test
    public void testSelection() {
        track.setSelectedSamples(Arrays.asList("D66", "CC-124", "2137"));
        assertEquals(Arrays.asList("2137", "CC-124", "D66"), visibleSamples());

        track.setSelectedSamples(null);
        assertEquals(15, track.sampleCount());
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
        assertEquals(false, json.has("selectedSamples"));
    }

    /**
     * Filtering by ID and by attribute are mutually exclusive -- choosing one turns the other off.
     */
    @Test
    public void testIdAndAttributeFiltersAreExclusive() {
        AttributeManager attributeManager = AttributeManager.getInstance();
        attributeManager.addAttribute("CC-124", "strain", "lab");
        attributeManager.addAttribute("CC-125", "strain", "lab");
        attributeManager.addAttribute("D66", "strain", "wild");
        SampleFilter labFilter = new SampleFilter(true,
                List.of(new FilterElement("strain", FilterElement.Operator.EQUAL, "lab")));

        track.setSelectedSamples(Arrays.asList("CC-124", "D66"));
        track.setSampleFilter(labFilter);
        assertNull(track.getSelectedSamples());
        assertEquals(Arrays.asList("CC-124", "CC-125"), visibleSamples());

        track.setSelectedSamples(Arrays.asList("CC-124", "D66"));
        assertNull(track.getSampleFilter());
        assertEquals(Arrays.asList("CC-124", "D66"), visibleSamples());

        // Clearing one filter does not restore the other
        track.setSelectedSamples(null);
        assertNull(track.getSampleFilter());
        assertEquals(15, track.sampleCount());
    }

    @Test
    public void testSessionRoundTrip() throws Exception {
        track.setSelectedSamples(Arrays.asList("CC-125", "cw15"));

        JSONObject json = new JSONObject();
        track.marshalJSON(json);

        VariantTrack restored = loadVariantTrack();
        restored.unmarshalJSON(json);
        assertEquals(Arrays.asList("CC-125", "cw15"), restored.getSelectedSamples());
        assertEquals(Arrays.asList("CC-125", "cw15"), visibleSamples(restored));
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
