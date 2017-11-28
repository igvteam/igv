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

import org.broad.igv.ui.javafx.JavaFXUIUtilities;

import javafx.scene.control.ScrollPane;
import javafx.scene.control.ScrollPane.ScrollBarPolicy;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;

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
        scrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
        scrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);

        JavaFXUIUtilities.bindWidthToContainer(mainContentPane, scrollPane);
        
        JavaFXUIUtilities.bindWidthToProperty(namePane, mainContentPane.namePaneWidthProperty());
        JavaFXUIUtilities.bindHeightToContainer(this, namePane);

        JavaFXUIUtilities.bindWidthToProperty(attributePane, mainContentPane.attributePaneWidthProperty());
        JavaFXUIUtilities.bindHeightToContainer(this, attributePane);

        // The contentContainer should take the rest of the space.  That is:
        // total width - (name pane width + attr pane width + (2 * insets) + scrollbar width)
        contentContainer.prefWidthProperty().bind(this.prefWidthProperty()
                .subtract(mainContentPane.namePaneWidthProperty()
                        .add(mainContentPane.attributePaneWidthProperty()).add(2 * INSET_SPACING + 30)));
        JavaFXUIUtilities.bindHeightToContainer(this, contentContainer);
        
        getChildren().add(namePane);
        getChildren().add(attributePane);
        getChildren().add(contentContainer);

        JavaFXUIUtilities.bindWidthToContainer(mainContentPane, this);

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
