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

// Intended as the rough equivalent of the IGVPanel class of the Swing UI, subclassed for handling the Header.  Work in progress.
public class HeaderRow extends IGVRow<NameHeaderPane, AttributeHeaderPane, HeaderPaneContainer> {

    public HeaderRow(MainContentPane mainContentPane) {
        init(mainContentPane, new NameHeaderPane(), new AttributeHeaderPane(), new HeaderPaneContainer());

        // Temporarily hard-coding the height since we have no header contents yet.
        int headerHeight = 100;
        getNamePane().prefHeightProperty().set(headerHeight);
        getAttributePane().prefHeightProperty().set(headerHeight);
        getContentContainer().prefHeightProperty().set(headerHeight);

        // TODO: move to CSS file
        getNamePane().setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(0, 0, 0)");
        getAttributePane().setStyle("-fx-border-style: solid; -fx-border-insets: 2; -fx-border-color: rgb(0, 0, 0)");
    }
}
