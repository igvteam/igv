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
import javafx.geometry.Orientation;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.ScrollPane.ScrollBarPolicy;
import javafx.scene.control.SplitPane;
import javafx.scene.layout.BorderPane;
import org.apache.log4j.Logger;
import org.broad.igv.prefs.PreferencesManager;

import static org.broad.igv.prefs.Constants.NAME_PANEL_WIDTH;
import static org.broad.igv.prefs.Constants.SHOW_SINGLE_TRACK_PANE_KEY;

// Intended as the rough equivalent of the MainPanel class of the Swing UI.  Work in progress.
public class MainContentPane extends BorderPane {
    private static Logger log = Logger.getLogger(MainContentPane.class);

    // Probably most/all components belong here
    private TrackScrollPane featureTrackScrollPane;
    private SplitPane centerSplitPane;

    // This should be attached directly to the PreferencesManager somehow
    private DoubleProperty namePaneWidthProp = new SimpleDoubleProperty(
            PreferencesManager.getPreferences().getAsFloat(NAME_PANEL_WIDTH));
    private DoubleProperty namePaneHiddenProp = new SimpleDoubleProperty(0);
    private DoubleProperty attributePaneWidthProp = new SimpleDoubleProperty(20);

    public MainContentPane() {
        HeaderRow headerRowContainer = new HeaderRow(this);
        ScrollPane headerScrollPane = new ScrollPane(headerRowContainer);
        this.setTop(headerScrollPane);

        // Mimic the SB policy of the Swing UI
        headerScrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
        headerScrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);

        centerSplitPane = new SplitPane();
        centerSplitPane.setOrientation(Orientation.VERTICAL);

        this.setCenter(centerSplitPane);


        TrackRow dataTrackRowContainer = new TrackRow(this);
        TrackScrollPane dataTrackScrollPane = new TrackScrollPane(dataTrackRowContainer);
        dataTrackScrollPane.setFitToHeight(true);
        dataTrackScrollPane.setFitToWidth(true);

        // Temporary, to show pane location
        dataTrackRowContainer.getNamePane().setStyle("-fx-background-color: green");

        // Mimic the SB policy of the Swing UI
        dataTrackScrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
        dataTrackScrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);

        centerSplitPane.getItems().add(dataTrackScrollPane);

        TrackRow featureTrackRowContainer = null;
        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            featureTrackRowContainer = new TrackRow(this);
            featureTrackScrollPane = new TrackScrollPane(featureTrackRowContainer);
            featureTrackScrollPane.setFitToHeight(true);
            featureTrackScrollPane.setFitToWidth(true);

            // Temporary, to show pane location
            featureTrackRowContainer.getNamePane().setStyle("-fx-background-color: yellow");

            // Mimic the SB policy of the Swing UI
            featureTrackScrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
            featureTrackScrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);

            centerSplitPane.getItems().add(featureTrackScrollPane);

            centerSplitPane.setDividerPositions(0.8);
        }
    }

    public void addInitialTrackRows() {
        // Guard to prevent repeat call


        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {

            centerSplitPane.setDividerPositions(0.8);
        }
    }

    public void setDividerLocations(double[] fractions) {
        centerSplitPane.setDividerPositions(fractions);
    }

    public double[] getDividerLocations() {
        return centerSplitPane.getDividerPositions();
    }

    public DoubleProperty getNamePaneWidthProp() {
        return namePaneWidthProp;
    }

    public DoubleProperty getAttributePaneWidthProp() {
        return attributePaneWidthProp;
    }
}
