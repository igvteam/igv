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

import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.Priority;

// Intended as the rough equivalent of the IGVPanel class of the Swing UI.  Work in progress.
public class HeaderRow extends HBox implements IGVRow {
    // Back-pointer to parent for coordinated management of rows
    private MainContentPane mainContentPane;

    // Do these need parent classes / interfaces above them? Perhaps all
    // IGVRowComponents (instead of TrackRowComponents)
    private NameHeaderPane nameHeaderPane = new NameHeaderPane();
    private AttributeHeaderPane attributeHeaderPane = new AttributeHeaderPane();
    private HeaderPaneContainer headerPaneContainer = new HeaderPaneContainer();

    public HeaderRow(MainContentPane mainContentPane) {
        super(5);
        this.mainContentPane = mainContentPane;
        getChildren().add(nameHeaderPane);
        getChildren().add(attributeHeaderPane);
        getChildren().add(headerPaneContainer);

        // Temporarily hard-coding the height since we have no header contents yet.
        int headerHeight = 100;
        nameHeaderPane.prefHeightProperty().set(headerHeight);
        nameHeaderPane.prefWidthProperty().bind(mainContentPane.getNamePaneWidthProp());
        attributeHeaderPane.prefWidthProperty().bind(mainContentPane.getAttributePaneWidthProp());

        // Set the headerPaneContainer to fill out the content width.
        HBox.setHgrow(headerPaneContainer, Priority.ALWAYS);

        // Temporary, to show pane location
        headerPaneContainer.setStyle("-fx-background-color: red");
        nameHeaderPane.setStyle("-fx-background-color: red");
        attributeHeaderPane.setStyle("-fx-background-color: blue");

        prefWidthProperty().bind(mainContentPane.prefWidthProperty());
        minWidthProperty().bind(mainContentPane.minWidthProperty());
        maxWidthProperty().bind(mainContentPane.maxWidthProperty());
    }

    @Override
    public MainContentPane getMainContentPane() {
        return mainContentPane;
    }

    @Override
    public Pane getNamePane() {
        return nameHeaderPane;
    }

    @Override
    public Pane getAttributePane() {
        return attributeHeaderPane;
    }

    @Override
    public Pane getContentContainer() {
        return headerPaneContainer;
    }

    public NameHeaderPane getNameHeaderPane() {
        return nameHeaderPane;
    }

    public AttributeHeaderPane getAttributeHeaderPane() {
        return attributeHeaderPane;
    }

    public HeaderPaneContainer getHeaderPaneContainer() {
        return headerPaneContainer;
    }
}
