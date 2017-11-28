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
package org.broad.igv.ui.javafx;

import javafx.beans.property.DoubleProperty;
import javafx.scene.layout.Region;

public class JavaFXUIUtilities {
    public static final void bindWidthToContainer(Region container, Region child) {
        child.prefWidthProperty().bind(container.prefWidthProperty());
        child.minWidthProperty().bind(container.minWidthProperty());
        child.maxWidthProperty().bind(container.maxWidthProperty());
    }

    public static final void bindWidthToProperty(Region component, DoubleProperty widthProperty) {
        component.prefWidthProperty().bind(widthProperty);
        component.minWidthProperty().bind(widthProperty);
        component.maxWidthProperty().bind(widthProperty);
    }

    public static final void bindHeightToContainer(Region container, Region child) {
        child.prefHeightProperty().bind(container.prefHeightProperty());
        child.minHeightProperty().bind(container.minHeightProperty());
        child.maxHeightProperty().bind(container.maxHeightProperty());
    }
    
    public static final void bindDimensionsToContainer(Region container, Region child) {
        bindWidthToContainer(container, child);
        bindHeightToContainer(container, child);
    }
}
