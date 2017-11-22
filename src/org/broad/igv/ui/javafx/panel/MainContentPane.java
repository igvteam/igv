/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broad.igv.ui.javafx.panel;

import javafx.application.Platform;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.control.SplitPane;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.CornerRadii;
import javafx.scene.paint.Color;
import org.apache.log4j.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.javafx.IGVToolBarManager;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;

import static org.broad.igv.prefs.Constants.*;

// Intended as the rough equivalent of the MainPanel class of the Swing UI.  Work in progress.
public class MainContentPane extends BorderPane {
    private static Logger log = Logger.getLogger(MainContentPane.class);
    
    // Probably most/all components should be instance vars.  Will migrate as need arises.
    private HeaderRow headerRow = null;
    private TrackRow featureTrackRow = null;
    private TrackScrollPane featureTrackScrollPane = null;
    private TrackRow dataTrackRow = null;
    private SplitPane centerSplitPane;

    private DoubleProperty namePaneWidthProp = new SimpleDoubleProperty(
            PreferencesManager.getPreferences().getAsFloat(NAME_PANEL_WIDTH));
    // TODO: find the proper sizing for the attribute panels
    private DoubleProperty attributePaneWidthProp = new SimpleDoubleProperty(20);
    private DoubleProperty dataPaneWidthProp = new SimpleDoubleProperty();
    private final Map<String, TrackRow> trackRowByName = new HashMap<String, TrackRow>();

    private IGVToolBarManager igvToolBarManager;
    
    public MainContentPane() {
    }

    public IGVToolBarManager getIgvToolBarManager() {
        return igvToolBarManager;
    }

    public void setIgvToolBarManager(IGVToolBarManager igvToolBarManager) {
        this.igvToolBarManager = igvToolBarManager;
    }

    // Used by the callback method of IGVStageBuilder to finish setting up this UI, after Stage init is done.
    public void initializeUI() {
        // Guard to prevent repeat call
        if (headerRow != null) {
            return;
        }

        dataPaneWidthProp.bind(this.prefWidthProperty().subtract(namePaneWidthProp).subtract(attributePaneWidthProp));

        headerRow = new HeaderRow(this);
        this.setTop(headerRow.getScrollPane());
        
        centerSplitPane = new SplitPane();
        centerSplitPane.setOrientation(Orientation.VERTICAL);

        this.setCenter(centerSplitPane);

        Color bgColor = PreferencesManager.getPreferences().getAsJavaFxColor(BACKGROUND_COLOR);
        Background background = new Background(new BackgroundFill(bgColor, CornerRadii.EMPTY, Insets.EMPTY));

        this.backgroundProperty().set(background);

        dataTrackRow = addTrackRow(IGV.DATA_PANEL_NAME);

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            featureTrackRow = addTrackRow(IGV.FEATURE_PANEL_NAME);
            featureTrackScrollPane = featureTrackRow.getScrollPane();

            centerSplitPane.setDividerPositions(0.9);
        }
    }


    // The following should only be called within Platform.runLater() or otherwise on the Application thread
    public TrackRow addTrackRow(String name) {
        TrackRow trackRow = new TrackRow(name, this);
        TrackScrollPane trackScrollPane = new TrackScrollPane(trackRow);
        trackRowByName.put(name, trackRow);

        int featurePaneIdx = -1;
        if (featureTrackScrollPane != null) {
            featurePaneIdx = centerSplitPane.getItems().indexOf(featureTrackScrollPane);
        }

        if (featurePaneIdx > 0) {
            centerSplitPane.getItems().add(featurePaneIdx, trackScrollPane);
        } else {
            centerSplitPane.getItems().add(trackScrollPane);
        }

        // Probably need to deal with centerSplitPane divider positions here.

        return trackRow;
    }
    
    
    public DoubleProperty namePaneWidthProperty() {
        return namePaneWidthProp;
    }

    public DoubleProperty attributePaneWidthProperty() {
        return attributePaneWidthProp;
    }

    public DoubleProperty dataPaneWidthProperty() {
        return dataPaneWidthProp;
    }

    public void hideNamePanel() {
        namePaneWidthProp.set(0);
    }

    public void showNamePanel() {
        namePaneWidthProp.set(PreferencesManager.getPreferences().getAsFloat(NAME_PANEL_WIDTH));
    }

    public boolean isNamePanelHidden() {
        return namePaneWidthProp.get() <= 0;
    }

    public Collection<TrackRow> getAllTrackRows() {
        return trackRowByName.values();
    }

    public HeaderRow getHeaderRow() {
        return headerRow;
    }
    
    public void resetTrackRows() {
        for (TrackRow trackRow : getAllTrackRows()) {
            trackRow.clearTracks();
            if (trackRow != featureTrackRow && trackRow != dataTrackRow) {
                centerSplitPane.getItems().remove(trackRow.getScrollPane());
            }
        }
        trackRowByName.clear();
        if (featureTrackRow != null) {
            trackRowByName.put(IGV.FEATURE_PANEL_NAME, featureTrackRow);
        }
        if (dataTrackRow != null) {
            trackRowByName.put(IGV.DATA_PANEL_NAME, dataTrackRow);
            dataTrackRow.reset();
        }
    }
    
    public TrackRow getTrackRow(String name) {
        TrackRow row = trackRowByName.get(name);
        if (row != null) {
            return row;
        }

        // If we get this far this is a new row
        FutureTask<TrackRow> trackRowCreator = new FutureTask<TrackRow>(new Callable<TrackRow>() {
            @Override
            public TrackRow call() throws Exception {
                return addTrackRow(name);
            }
        });
        
        Platform.runLater(trackRowCreator);
        try {
            return trackRowCreator.get();
        }
        catch (ExecutionException | InterruptedException e) {
            // TODO: better error handling.  Prob need equivalent of MessageUtils.showMessage() 
            log.error(e);
            throw new RuntimeException(e);
        }
    }
    
    // Not yet used; anticipated...
    public void setDividerLocations(double[] fractions) {
        centerSplitPane.setDividerPositions(fractions);
    }

    public double[] getDividerLocations() {
        return centerSplitPane.getDividerPositions();
    }
}
