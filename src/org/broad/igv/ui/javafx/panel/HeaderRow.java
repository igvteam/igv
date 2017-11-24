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

import javafx.scene.control.ScrollPane;

// Intended as the rough equivalent of the IGVPanel class of the Swing UI, subclassed for handling the Header.  Work in progress.
public class HeaderRow extends IGVRow<NameHeaderPane, AttributeHeaderPane, HeaderPaneContainer, ScrollPane> {

    public HeaderRow(MainContentPane mainContentPane) {
        // TODO: determine correct header height setting
        double headerHeight = 140.0;
        this.prefHeightProperty().set(headerHeight);
        ScrollPane scrollPane = new ScrollPane(this);
        NameHeaderPane namePane = new NameHeaderPane();
        AttributeHeaderPane attributePane = new AttributeHeaderPane();
        HeaderPaneContainer headerPaneContainer = new HeaderPaneContainer();
        init(mainContentPane, namePane, attributePane, headerPaneContainer, scrollPane);

        // TODO: move to CSS file
        scrollPane.setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(102, 102, 102)");
        getNamePane().setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(0, 0, 0)");
        getAttributePane().setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(0, 0, 0)");

        // Wait to create the headerPaneContainer content until after layout & sizing is complete.
        headerPaneContainer.createHeaderPanes();
    }
}
