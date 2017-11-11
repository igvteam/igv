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

import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.Node;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.ScrollPane.ScrollBarPolicy;
import javafx.scene.control.SplitPane;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.CornerRadii;
import javafx.scene.paint.Color;
import org.apache.commons.lang.StringUtils;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;

import java.util.ArrayList;
import java.util.List;

import static org.broad.igv.prefs.Constants.*;

// Intended as the rough equivalent of the MainPanel class of the Swing UI.  Work in progress.
public class MainContentPane extends BorderPane {

    // Probably most/all components should be instance vars.  Will migrate as need arises.
    private TrackRow featureTrackRowContainer = null;
    private TrackScrollPane featureTrackScrollPane = null;
    private SplitPane centerSplitPane;

    private DoubleProperty namePaneWidthProp = new SimpleDoubleProperty(
            PreferencesManager.getPreferences().getAsFloat(NAME_PANEL_WIDTH));
    private DoubleProperty attributePaneWidthProp = new SimpleDoubleProperty(20);


    public MainContentPane() {
    }

    // Used by the callback method of IGVStageBuilder to finish setting up this UI, after Stage init is done.
    public void initializeUI() {

        // Need guard to prevent repeat call

        HeaderRow headerRow = new HeaderRow(this);
        ScrollPane headerScrollPane = new ScrollPane(headerRow);
        headerScrollPane.setFitToHeight(true);
        headerScrollPane.setFitToWidth(true);

        // TODO: move to CSS file
        headerScrollPane.setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(102, 102, 102)");
        this.setTop(headerScrollPane);

        // Mimic the SB policy of the Swing UI
        headerScrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
        headerScrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);

        centerSplitPane = new SplitPane();
        centerSplitPane.setOrientation(Orientation.VERTICAL);

        this.setCenter(centerSplitPane);

        Color bgColor = PreferencesManager.getPreferences().getAsJavaFxColor(BACKGROUND_COLOR);
        Background background = new Background(new BackgroundFill(bgColor, CornerRadii.EMPTY, Insets.EMPTY));

        this.backgroundProperty().set(background);

        TrackRow dataTrackRow = new TrackRow(IGV.DATA_PANEL_NAME, this);
        TrackScrollPane dataTrackScrollPane = new TrackScrollPane(dataTrackRow);

        centerSplitPane.getItems().add(dataTrackScrollPane);

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            featureTrackRowContainer = new TrackRow(IGV.FEATURE_PANEL_NAME, this);
            featureTrackScrollPane = new TrackScrollPane(featureTrackRowContainer);

            centerSplitPane.getItems().add(featureTrackScrollPane);

            centerSplitPane.setDividerPositions(0.9);
        }
    }

    // The following should be called within Platform.runLater()
    public TrackScrollPane addTrackRow(String name) {
        TrackRow dataTrackRow = new TrackRow(name, this);
        TrackScrollPane dataTrackScrollPane = new TrackScrollPane(dataTrackRow);

        int featurePaneIdx = -1;
        if (featureTrackScrollPane != null) {
            featurePaneIdx = centerSplitPane.getItems().indexOf(featureTrackScrollPane);
        }

        if (featurePaneIdx > 0) {
            centerSplitPane.getItems().add(featurePaneIdx, dataTrackScrollPane);
        } else {
            centerSplitPane.getItems().add(dataTrackScrollPane);
        }

        return dataTrackScrollPane;
    }
    
    
    public DoubleProperty namePaneWidthProperty() {
        return namePaneWidthProp;
    }

    public DoubleProperty attributePaneWidthProperty() {
        return attributePaneWidthProp;
    }

    public void hideNamePanel() {
        namePaneWidthProp.set(0);
    }

    public void showNamePanel() {
        namePaneWidthProp.set(PreferencesManager.getPreferences().getAsFloat(NAME_PANEL_WIDTH));
    }

    public List<TrackRow> getAllTrackRows() {
        // Just porting Swing code for now; seems clumsy to build this every time and query via instanceof
        // TODO: consider keeping an instance var list that we manage during addTrackRow() etc.
        // May need to synchronize access.  Maybe should be a Map to help out with getTrackRow(name)
        // Or just keep a List of Tracks?  Have to find what exactly is needed.
        List<TrackRow> trackRows = new ArrayList<TrackRow>();
        for (Node child : centerSplitPane.getChildrenUnmodifiable()) {
            if (child instanceof TrackScrollPane) {
                trackRows.add(((TrackScrollPane) child).getTrackRow());
            }
        }

        return trackRows;
    }

    public TrackRow getTrackRow(String name) {
        List<TrackRow> trackRows = getAllTrackRows();
        for (TrackRow row : trackRows) {
            if (StringUtils.equals(name, row.getName())) {
                return row;
            }
        }

        // If we get this far this is a new row
        TrackScrollPane trackScrollPane = addTrackRow(name);
        return trackScrollPane.getTrackRow();
    }
    
    public boolean isNamePanelHidden() {
        return namePaneWidthProp.get() <= 0;
    }

    // Not yet used; anticipated...
    public void setDividerLocations(double[] fractions) {
        centerSplitPane.setDividerPositions(fractions);
    }

    public double[] getDividerLocations() {
        return centerSplitPane.getDividerPositions();
    }
}
