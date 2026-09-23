package org.igv.track;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.BasicFeature;
import org.igv.track.Track.DisplayMode;
import org.igv.ui.panel.TrackPanelScrollPane;
import org.json.JSONObject;
import org.junit.Test;

import java.util.Collections;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests the interaction between display mode, row height, and track height for tracks with rows.
 */
public class DisplayModeRowHeightTest extends AbstractHeadlessTest {

    private FeatureTrack newFeatureTrack() {
        FeatureCollectionSource<BasicFeature> source = new FeatureCollectionSource<>(Collections.emptyList(), genome);
        return new FeatureTrack("test", "test", source);
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
