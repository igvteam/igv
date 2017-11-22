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
import javafx.scene.control.ScrollPane.ScrollBarPolicy;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.Priority;

// Intended as the rough equivalent of the IGVPanel class of the Swing UI.  Work in progress.
// Note: N, A, C might need to become more specific types later (e.g. RowComponent).  Pane is general enough for now.
public class IGVRow<N extends Pane, A extends Pane, C extends Pane, S extends ScrollPane> extends HBox {

    private static final double INSET_SPACING = 5;

    private MainContentPane mainContentPane;
    private N namePane;
    private A attributePane;
    private C contentContainer;
    private S scrollPane;

    protected IGVRow() {
        super(INSET_SPACING);
    }

    protected void init(MainContentPane mainContentPane, N namePane, A attributePane, C contentContainer, S scrollPane) {

        this.mainContentPane = mainContentPane;
        this.namePane = namePane;
        this.attributePane = attributePane;
        this.contentContainer = contentContainer;
        this.scrollPane = scrollPane;

        scrollPane.setFitToHeight(true);
        scrollPane.setFitToWidth(true);

        // Mimic the SB policy of the Swing UI
        scrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
        scrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);
        
        namePane.prefWidthProperty().bind(mainContentPane.namePaneWidthProperty());
        namePane.minWidthProperty().bind(mainContentPane.namePaneWidthProperty());
        namePane.maxWidthProperty().bind(mainContentPane.namePaneWidthProperty());
        namePane.prefHeightProperty().bind(this.prefHeightProperty());
        namePane.minHeightProperty().bind(this.minHeightProperty());
        namePane.maxHeightProperty().bind(this.maxHeightProperty());

        attributePane.prefWidthProperty().bind(mainContentPane.attributePaneWidthProperty());
        attributePane.minWidthProperty().bind(mainContentPane.attributePaneWidthProperty());
        attributePane.maxWidthProperty().bind(mainContentPane.attributePaneWidthProperty());
        attributePane.prefHeightProperty().bind(this.prefHeightProperty());
        attributePane.minHeightProperty().bind(this.minHeightProperty());
        attributePane.maxHeightProperty().bind(this.maxHeightProperty());

        contentContainer.prefWidthProperty().bind(this.prefWidthProperty()
                .subtract(mainContentPane.namePaneWidthProperty())
                .subtract(mainContentPane.attributePaneWidthProperty()));
        
        getChildren().add(namePane);
        getChildren().add(attributePane);
        getChildren().add(contentContainer);

        // Set the contentContainer to fill out the content width.
        HBox.setHgrow(contentContainer, Priority.ALWAYS);

        prefWidthProperty().bind(mainContentPane.prefWidthProperty());
        minWidthProperty().bind(mainContentPane.minWidthProperty());
        maxWidthProperty().bind(mainContentPane.maxWidthProperty());

        backgroundProperty().bind(mainContentPane.backgroundProperty());
        namePane.backgroundProperty().bind(mainContentPane.backgroundProperty());
        attributePane.backgroundProperty().bind(mainContentPane.backgroundProperty());
        contentContainer.backgroundProperty().bind(mainContentPane.backgroundProperty());
    }

    public MainContentPane getMainContentPane() {
        return mainContentPane;
    }

    public N getNamePane() {
        return namePane;
    }

    public A getAttributePane() {
        return attributePane;
    }

    public C getContentContainer() {
        return contentContainer;
    }

    public S getScrollPane() {
        return scrollPane;
    }
}
