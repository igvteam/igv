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
        computeFrameBounds();
        for (ReferenceFrame f : frames) {
            if (f.isVisible()) {
                log.info("creating HeaderPane for " + f.getChrName());
                HeaderPane headerPane = new HeaderPane(f);
                headerPanes.add(headerPane);
                // TODO: Need to account for multiple frames in width.  The following is wrong.
                // We need to split the width among all HPs.  Prob extract from the RefFrame?
                headerPane.prefWidthProperty().bind(prefWidthProperty());
                headerPane.minWidthProperty().bind(minWidthProperty());
                headerPane.maxWidthProperty().bind(maxWidthProperty());
                headerPane.backgroundProperty().bind(backgroundProperty());
                headerPane.prefHeightProperty().bind(prefHeightProperty());
                headerPane.minHeightProperty().bind(minHeightProperty());
                headerPane.maxHeightProperty().bind(maxHeightProperty());
                contentPane.getChildren().add(headerPane);
            }
        }
        contentPane.prefWidthProperty().bind(prefWidthProperty());
        
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

        this.prefWidthProperty().addListener((observable, oldValue, newValue) -> computeFrameBounds());
    }
    
    private void computeFrameBounds() {
        List<ReferenceFrame> frames = FrameManager.getFrames();
        Double width = prefWidthProperty().get();
        if (frames.size() == 1) {
            // TODO: convert bounds to double or DoubleProperty for JavaFX?
            frames.get(0).setBounds(0, width.intValue());
        }
        else {
            // Not yet...
        }

        // Here's a thought on splitting the width equally:
        // Create a contentPrefWidthProperty to bind to each of the headerPane's prefWidthProps.
        // It would be bound as: prop.bind(this.prefWidthProperty().divide(frames.size());
        // For the usual case of size == 1, just bind directly for better perf.
        // Then we would just need to set the pixelX value for each RefFrame.
        // Perhaps the RefFrame could just work based on the origin?  Seems like we can get the
        // pixelX from the HeaderPane (again as a property).
        // 
        // In any case, we'll go this route for now to prove it out.
    }
}
