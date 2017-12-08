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

import javafx.event.EventHandler;
import javafx.scene.control.Label;
import javafx.scene.control.Tooltip;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.Pane;
import org.apache.log4j.Logger;
import org.broad.igv.ui.javafx.IGVBackendPlaceholder;
import org.broad.igv.ui.javafx.JavaFXUIUtilities;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

// Intended as the rough equivalent of the HeaderPanel class of the Swing UI.  Work in progress.
// TODO: Need to add equivalents of RegionOfInterestPanel, GeneListPanel, etc.
// TODO: DnD handling
public class HeaderPane extends BorderPane {
    private static Logger log = Logger.getLogger(HeaderPane.class);
    
    private ReferenceFrame frame;
    private BorderPane geneListPane;

    private CytobandPane cytobandPane;
    private RulerPane rulerPane;

    public HeaderPane(ReferenceFrame frame) {

        this.frame = frame;
        JavaFXUIUtilities.bindWidthToProperty(this, frame.displayWidthProperty());

        if (FrameManager.isGeneListMode()) {
            Label label = new Label(this.frame.getName());
            label.setStyle("-fx-alignment: center; -fx-text-fill: blue;");
            label.prefWidthProperty().bind(this.prefWidthProperty());
            label.setTooltip(new Tooltip(frame.getName()));
            this.geneListPane = new BorderPane();
            cytobandPane = new CytobandPane(frame, false);
            geneListPane.setCenter(cytobandPane);
            geneListPane.setBottom(label);
            this.setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: gray; -fx-background-color: yellow;");

            JavaFXUIUtilities.bindWidthToContainer(geneListPane, label);
            cytobandPane.prefWidthProperty().bind(geneListPane.prefWidthProperty().subtract(4));
            cytobandPane.minWidthProperty().bind(geneListPane.minWidthProperty().subtract(4));
            cytobandPane.maxWidthProperty().bind(geneListPane.maxWidthProperty().subtract(4));

            geneListPane.prefWidthProperty().bind(this.prefWidthProperty().subtract(4));
            geneListPane.minWidthProperty().bind(this.minWidthProperty().subtract(4));
            geneListPane.maxWidthProperty().bind(this.maxWidthProperty().subtract(4));

            this.setOnMouseClicked(geneListModeOnClickHandler);
            cytobandPane.setOnMouseClicked(geneListModeOnClickHandler);
            label.setOnMouseClicked(geneListModeOnClickHandler);
            geneListPane.setOnMouseClicked(geneListModeOnClickHandler);
            
            this.setCenter(geneListPane);
        } else {
            BorderPane pane = new BorderPane();
            cytobandPane = new CytobandPane(frame);
            rulerPane = new RulerPane(frame);
            pane.setTop(cytobandPane);
            pane.setCenter(rulerPane);

            JavaFXUIUtilities.bindWidthToContainer(this, rulerPane);
            rulerPane.backgroundProperty().bind(backgroundProperty());
            JavaFXUIUtilities.bindWidthToContainer(this, cytobandPane);

            Pane roiDummy = new Pane();
            roiDummy.setMinSize(15.0, 15.0);
            pane.setBottom(roiDummy);
            this.setCenter(pane);
        }

        cytobandPane.backgroundProperty().bind(backgroundProperty());
    }

    private EventHandler<MouseEvent> geneListModeOnClickHandler = (event) -> {
        if (event.getClickCount() > 1) {
            FrameManager.setToDefaultFrame(frame.getName());
            IGVBackendPlaceholder.resetContent();
            event.consume();
        }
    };
}
