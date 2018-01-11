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
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.igv.lists.GeneList;
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
    private HBox contentPane = new HBox();
    
    public HeaderPaneContainer() {
        JavaFXUIUtilities.bindWidthToContainer(this, contentPane);
        contentPane.prefHeightProperty().bind(prefHeightProperty());
        setCenter(contentPane);
    }

    public void createHeaderPanes() {
        headerPanes.clear();
        contentPane.getChildren().clear();

        List<ReferenceFrame> frames = FrameManager.getFrames();
        for (ReferenceFrame f : frames) {
            if (f.isVisible()) {
                HeaderPane headerPane = new HeaderPane(f);
                headerPanes.add(headerPane);
                headerPane.backgroundProperty().bind(backgroundProperty());
                JavaFXUIUtilities.bindHeightToContainer(this, headerPane);
                contentPane.getChildren().add(headerPane);
            }
        }

        if (FrameManager.isGeneListMode()) {
            GeneList gl = IGVBackendPlaceholder.getCurrentGeneList();
            String name = gl.getDisplayName();
            if (StringUtils.isNotBlank(name)) {
                Label label = new Label(name);
                label.getStyleClass().add("geneListHeaderPaneContainerLabel");
                contentPane.prefHeightProperty().bind(prefHeightProperty().subtract(label.heightProperty()));
                setTop(label);
            }
        }
    }

    public DoubleProperty frameSpacingProperty() {
        return contentPane.spacingProperty();
    }
}
