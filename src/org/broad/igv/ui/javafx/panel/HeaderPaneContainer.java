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
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

// Intended as the rough equivalent of the HeaderPanelContainer class of the Swing UI.  Work in progress.
public class HeaderPaneContainer extends BorderPane {

    public HeaderPaneContainer() {
        //setStyle("-fx-background-color: purple; -fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(0, 0, 0); -fx-backgound-color: red");
        createHeaderPanes();
    }

    public void createHeaderPanes() {
        getChildren().removeAll();

        HBox contentPane = new HBox(6);
        contentPane.setStyle("-fx-background-color: purple");
        for (ReferenceFrame f : FrameManager.getFrames()) {
            if (f.isVisible()) {
                HeaderPane headerPane = new HeaderPane(f);
//                headerPane.prefHeightProperty().bind(prefHeightProperty());
//                headerPane.prefWidthProperty().bind(prefWidthProperty());
//                headerPane.backgroundProperty().bind(backgroundProperty());
                headerPane.setStyle("-fx-background-color: green");
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
            // Display label for testing
            String name = "Test";
            Label label = new Label(name);
            label.setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: lightgray");
            // Normally will not subtract here
            contentPane.prefHeightProperty().bind(prefHeightProperty().subtract(label.heightProperty()));
            setTop(label);
        }

        setCenter(contentPane);
    }
}
