package org.igv.track;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.BasicFeature;
import org.igv.track.Track.DisplayMode;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.panel.TrackPanelScrollPane;
import org.json.JSONObject;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests the interaction between display mode, row height, and track height for tracks with rows.
 */
public class DisplayModeRowHeightTest extends AbstractHeadlessTest {

    private FeatureTrack newFeatureTrack() {
        return newFeatureTrack(Collections.emptyList());
    }

    private FeatureTrack newFeatureTrack(List<BasicFeature> features) {
        return new FeatureTrack("test", "test", new FeatureCollectionSource<>(features, genome));
    }

    /**
     * Eight overlapping features on chr1, so each packs into its own row.
     */
    private static List<BasicFeature> overlappingFeatures() {
        List<BasicFeature> features = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            features.add(new BasicFeature("chr1", 1000 + 10 * i, 1100));
        }
        return features;
    }

    private static ReferenceFrame frameOverFeatures() {
        ReferenceFrame frame = new ReferenceFrame("test");
        frame.setBounds(0, 1000);
        frame.jumpTo("chr1", 900, 1200);
        return frame;
    }

    @Test
    public void testDisplayModeSetsDefaultRowHeight() {
        FeatureTrack track = newFeatureTrack();
        assertTrue(track.hasRows());

        track.setDisplayMode(DisplayMode.SQUISHED);
        assertEquals(track.getDefaultSquishedRowHeight(), track.getRowHeight());

        track.setDisplayMode(DisplayMode.EXPANDED);
        assertEquals(track.getDefaultExpandedRowHeight(), track.getRowHeight());

        // An explicit row height overrides the default until the mode changes again
        track.setRowHeight(7);
        assertEquals(7, track.getRowHeight());
        track.setDisplayMode(DisplayMode.SQUISHED);
        assertEquals(track.getDefaultSquishedRowHeight(), track.getRowHeight());
    }

    @Test
    public void testModeChangeDoesNotChangeDisplayedTrackHeight() {
        FeatureTrack track = newFeatureTrack();
        track.setDisplayMode(DisplayMode.EXPANDED);

        // Not yet displayed -- the track remains auto-sized to its content
        track.setDisplayMode(DisplayMode.SQUISHED);
        assertEquals(track.getContentHeight(), track.getHeight());
        track.setDisplayMode(DisplayMode.EXPANDED);
        assertEquals(track.getContentHeight(), track.getHeight());

        // Once displayed, changing the mode changes the content height but not the track height
        track.setViewport(new TrackPanelScrollPane());
        int expandedHeight = track.getHeight();
        int expandedContentHeight = track.getContentHeight();
        track.setDisplayMode(DisplayMode.SQUISHED);
        assertTrue(track.getContentHeight() < expandedContentHeight);
        assertEquals(expandedHeight, track.getHeight());
        track.setDisplayMode(DisplayMode.EXPANDED);
        assertEquals(expandedHeight, track.getHeight());
    }

    @Test
    public void testModeChangeFitsHeightToContent() {
        FeatureTrack track = newFeatureTrack(overlappingFeatures());
        track.setDisplayMode(DisplayMode.COLLAPSED);
        track.load(frameOverFeatures());
        track.setViewport(new TrackPanelScrollPane());
        int collapsedHeight = track.getHeight();

        // Expand grows the track to the content height, but no further than the maximum
        track.setDisplayMode(DisplayMode.EXPANDED);
        track.adjustHeightForDisplayMode();
        assertTrue(track.getContentHeight() > track.getMaxExpandedHeight());
        assertEquals(track.getMaxExpandedHeight(), track.getHeight());

        // Squish shrinks the track to the content height
        track.setDisplayMode(DisplayMode.SQUISHED);
        track.adjustHeightForDisplayMode();
        assertTrue(track.getContentHeight() > collapsedHeight);
        assertEquals(track.getContentHeight(), track.getHeight());

        // Expand restores the height a squish or collapse shrank from, even above the type maximum
        track.setHeight(250);
        track.setDisplayMode(DisplayMode.COLLAPSED);
        track.adjustHeightForDisplayMode();
        assertEquals(collapsedHeight, track.getHeight());
        track.setDisplayMode(DisplayMode.EXPANDED);
        track.adjustHeightForDisplayMode();
        assertEquals(250, track.getHeight());

        // An explicit height is kept
        track.setExplicitHeight(50);
        track.setDisplayMode(DisplayMode.EXPANDED);
        track.adjustHeightForDisplayMode();
        assertEquals(50, track.getHeight());

        // Set Track Height 0 means auto-size, not an explicit height
        track.setExplicitHeight(0);
        track.setDisplayMode(DisplayMode.SQUISHED);
        track.adjustHeightForDisplayMode();
        assertEquals(track.getContentHeight(), track.getHeight());
    }

    @Test
    public void testFitDeferredUntilFeaturesLoad() {
        FeatureTrack track = newFeatureTrack(overlappingFeatures());
        track.setDisplayMode(DisplayMode.COLLAPSED);
        track.setViewport(new TrackPanelScrollPane());
        int collapsedHeight = track.getHeight();

        // Nothing loaded yet: the mode change pins the height and the fit waits for the load
        track.setDisplayMode(DisplayMode.EXPANDED);
        track.adjustHeightForDisplayMode();
        assertEquals(collapsedHeight, track.getHeight());

        track.load(frameOverFeatures());
        assertEquals(track.getMaxExpandedHeight(), track.getHeight());
    }

    @Test
    public void testCustomRowHeight() {
        FeatureTrack track = newFeatureTrack();
        track.setDisplayMode(DisplayMode.SQUISHED);

        // A custom row height leaves SQUISHED, so that SQUISHED can be selected again to restore its default
        track.setCustomRowHeight(40);
        assertEquals(DisplayMode.CUSTOM, track.getDisplayMode());
        assertEquals(40, track.getRowHeight());
        track.setDisplayMode(DisplayMode.SQUISHED);
        assertEquals(track.getDefaultSquishedRowHeight(), track.getRowHeight());

        // A displayed, auto-sized track resizes to the new row height rather than being pinned by the mode change
        track.setViewport(new TrackPanelScrollPane());
        track.setCustomRowHeight(60);
        assertEquals(track.getContentHeight(), track.getHeight());

        // COLLAPSED is a packing mode, not a row height, and is kept
        track.setDisplayMode(DisplayMode.COLLAPSED);
        track.setCustomRowHeight(40);
        assertEquals(DisplayMode.COLLAPSED, track.getDisplayMode());
    }

    @Test
    public void testCustomSessionRoundTrip() {
        FeatureTrack track = newFeatureTrack();
        track.setDisplayMode(DisplayMode.EXPANDED);
        track.setCustomRowHeight(40);

        // Written as EXPANDED for igv.js, which has no CUSTOM mode
        JSONObject json = new JSONObject();
        track.marshalJSON(json);
        assertEquals("EXPANDED", json.getString("displayMode"));

        FeatureTrack restored = newFeatureTrack();
        restored.unmarshalJSON(json);
        assertEquals(DisplayMode.CUSTOM, restored.getDisplayMode());
        assertEquals(40, restored.getRowHeight());

        // A default row height restores the plain mode
        track.setDisplayMode(DisplayMode.SQUISHED);
        json = new JSONObject();
        track.marshalJSON(json);
        restored = newFeatureTrack();
        restored.unmarshalJSON(json);
        assertEquals(DisplayMode.SQUISHED, restored.getDisplayMode());
    }
}
