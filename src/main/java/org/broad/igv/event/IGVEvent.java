package org.broad.igv.event;

import org.broad.igv.oauth.OAuthProvider;
import org.broad.igv.prefs.PreferencesChangeEvent;
import org.broad.igv.sam.InsertionSelectionEvent;
import org.broad.igv.ui.panel.FrameManager;

public sealed interface IGVEvent
        permits AlignmentTrackEvent, DataLoadedEvent, GenomeChangeEvent, GenomeResetEvent, RefreshEvent, StopEvent,
        TrackFilterEvent, TrackGroupEvent, ViewChange, OAuthProvider.AuthStateEvent, PreferencesChangeEvent,
        InsertionSelectionEvent, FrameManager.ChangeEvent {

}
