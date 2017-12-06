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
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.igv.lists.GeneList;
import org.broad.igv.session.Session;
import org.broad.igv.ui.javafx.IGVBackendPlaceholder;
import org.broad.igv.ui.javafx.JavaFXUIUtilities;
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
        int frameCount = frames.size();
        if (frameCount == 1) {
            HeaderPane headerPane = new HeaderPane(frames.get(0));
            headerPanes.add(headerPane);
            JavaFXUIUtilities.bindWidthToContainer(this, headerPane);
            headerPane.backgroundProperty().bind(backgroundProperty());
            JavaFXUIUtilities.bindHeightToContainer(this, headerPane);
            contentPane.getChildren().add(headerPane);
        } else {
            DoubleProperty widthForFrames = new SimpleDoubleProperty();
            widthForFrames.bind(this.prefWidthProperty().subtract(((frameCount - 1) * 6)).divide(frameCount));
            for (ReferenceFrame f : frames) {
                if (f.isVisible()) {
                    HeaderPane headerPane = new HeaderPane(f);
                    headerPanes.add(headerPane);
                    // TODO: Need to account for multiple frames in width.  The following is wrong.
                    // We need to split the width among all HPs.  Prob extract from the RefFrame?
//                    headerPane.setPrefWidth(f.getWidthInPixels());
//                    headerPane.setMinWidth(f.getWidthInPixels());
//                    headerPane.setMaxWidth(f.getWidthInPixels());
                    JavaFXUIUtilities.bindWidthToProperty(headerPane, widthForFrames);
                    headerPane.backgroundProperty().bind(backgroundProperty());
                    JavaFXUIUtilities.bindHeightToContainer(this, headerPane);
                    contentPane.getChildren().add(headerPane);
                }
            }
        }
        JavaFXUIUtilities.bindWidthToContainer(this, contentPane);

        contentPane.prefHeightProperty().bind(prefHeightProperty());
        if (FrameManager.isGeneListMode()) {
            GeneList gl = IGVBackendPlaceholder.getCurrentGeneList();
            String name = gl.getDisplayName();
            if (StringUtils.isNotBlank(name)) {
                Label label = new Label(name);
                label.setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: lightgray; -fx-alignment: center;");
                contentPane.prefHeightProperty().bind(prefHeightProperty().subtract(label.heightProperty()));
                setTop(label);
            }
        }

        setCenter(contentPane);

        this.prefWidthProperty().addListener((observable, oldValue, newValue) -> computeFrameBounds());
    }

    private int hgap = 5;

    private void computeFrameBounds() {
        List<ReferenceFrame> frames = FrameManager.getFrames();
        Double width = prefWidthProperty().get();
        if (frames.size() == 1) {
            // TODO: convert bounds to double or DoubleProperty for JavaFX?
            frames.get(0).setBounds(0, width.intValue());
        }
        else {

            float gap = Math.min(1, 20.0f / ((int) (1.5 * frames.size()))) * hgap;
            int x = 0;

            // Not dealing with Session for now.
            Session.GeneListMode mode = IGVBackendPlaceholder.getGeneListMode();
            double wc = mode == Session.GeneListMode.NORMAL ?
                    (width - (frames.size() - 1) * gap) / frames.size() : 20;

            for (int i = 0; i < frames.size(); i++) {
                ReferenceFrame frame = frames.get(i);
                int nextX = (int) ((i + 1) * (wc + gap));
                int w = nextX - x;
                frame.setBounds(x, w);
                x = nextX;
            }
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
