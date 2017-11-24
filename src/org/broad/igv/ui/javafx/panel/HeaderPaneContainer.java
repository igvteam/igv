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

import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import org.apache.log4j.Logger;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.ArrayList;
import java.util.List;

// Intended as the rough equivalent of the HeaderPanelContainer class of the Swing UI.  Work in progress.
public class HeaderPaneContainer extends BorderPane {
    private static Logger log = Logger.getLogger(HeaderPaneContainer.class);
    private List<HeaderPane> headerPanes = new ArrayList<HeaderPane>();
    
    public HeaderPaneContainer() {
    }

    public void createHeaderPanes() {
        getChildren().removeAll();
        headerPanes.clear();

        HBox contentPane = new HBox(6);

        List<ReferenceFrame> frames = FrameManager.getFrames();
        if (frames.size() == 1) {
            // TODO: convert bounds to double or DoubleProperty for JavaFX.
            double width = prefWidthProperty().get();
            frames.get(0).setBounds(0, (int) width);
        }

        for (ReferenceFrame f : frames) {
            if (f.isVisible()) {
                log.info("creating HeaderPane for " + f.getChrName());
                HeaderPane headerPane = new HeaderPane(f);
                headerPanes.add(headerPane);
                // Not correct; we need to split the width among all HPs
                // TODO: Need to account for multiple frames in width.  The following is wrong.
                headerPane.prefWidthProperty().bind(prefWidthProperty());
                headerPane.minWidthProperty().bind(minWidthProperty());
                headerPane.maxWidthProperty().bind(maxWidthProperty());
                headerPane.backgroundProperty().bind(backgroundProperty());
                headerPane.prefHeightProperty().bind(prefHeightProperty());
                headerPane.minHeightProperty().bind(minHeightProperty());
                headerPane.maxHeightProperty().bind(maxHeightProperty());

                log.info("HP HW: " + headerPane.getWidth() + ":" + headerPane.getHeight());
                log.info("HP pHW: " + headerPane.getPrefWidth() + ":" + headerPane.getPrefHeight());
                contentPane.getChildren().add(headerPane);
            }
        }
        contentPane.prefWidthProperty().bind(prefWidthProperty());

        // Here's a thought on splitting the width equally:
        // Create a contentPrefWidthProperty to bind to each of the headerPane's prefWidthProps.
        // It would be bound as: prop.bind(this.prefWidthProperty().divide(frames.size());
        // For the usual case of size == 1, just bind directly for better perf.
        
        if (FrameManager.isGeneListMode()) {
            GeneList gl = IGV.getInstance().getSession().getCurrentGeneList();
            String name = gl.getDisplayName();
            Label label = new Label(name);
            label.setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: lightgray");
            contentPane.prefHeightProperty().bind(prefHeightProperty().subtract(label.heightProperty()));
            setTop(label);
        } else {
            contentPane.prefHeightProperty().bind(prefHeightProperty());
        }

        setCenter(contentPane);
        log.info("HPC HW: " + getWidth() + ":" + getHeight());
        log.info("HPC pHW: " + getPrefWidth() + ":" + getPrefHeight());
        log.info("HPC content HW: " + contentPane.getWidth() + ":" + contentPane.getHeight());
        log.info("HPC content pHW: " + contentPane.getPrefWidth() + ":" + contentPane.getPrefHeight());
    }
}
