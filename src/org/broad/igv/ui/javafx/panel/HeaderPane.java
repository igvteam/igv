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
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

// Intended as the rough equivalent of the HeaderPanel class of the Swing UI.  Work in progress.
// TODO: Need to add equivalents of CytobandPanel, RulerPanel, RegionOfInterestPanel, GeneListPanel, etc.
// TODO: DnD handling
public class HeaderPane extends BorderPane {

    private ReferenceFrame frame;
    private BorderPane geneListPane;
    private Label label;

    private CytobandPane cytobandPane;

    public HeaderPane(ReferenceFrame frame) {
        // TODO: Set background, sizing, etc on all components

        this.frame = frame;

//        cytobandPane.prefHeightProperty().bind(prefHeightProperty());
//        cytobandPane.prefWidthProperty().bind(prefWidthProperty());
//        cytobandPane.backgroundProperty().bind(backgroundProperty());

        if (FrameManager.isGeneListMode()) {
            this.label = new Label(this.frame.getName());
            this.geneListPane = new BorderPane();
            setPrefSize(400, 100);
            cytobandPane = new CytobandPane(frame);
            geneListPane.setCenter(cytobandPane);
            geneListPane.setBottom(label);
            this.getChildren().add(geneListPane);
        } else {
            BorderPane pane = new BorderPane();
            cytobandPane = new CytobandPane(frame);
            pane.setTop(cytobandPane);
            this.getChildren().add(pane);
        }
        
    }

}
