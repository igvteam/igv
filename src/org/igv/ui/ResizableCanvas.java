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
package org.igv.ui;

import javafx.scene.canvas.Canvas;
import javafx.scene.layout.Pane;

// Canvas doesn't resize on its own, so we need to wrap it in a container to handle this on its behalf
// http://fxexperience.com/2014/05/resizable-grid-using-canvas/
// https://stackoverflow.com/questions/27808210/java-fx-splitpane
// Note that this is a simpler approach than presented above: just bind the H&W properties to the wrapper instead of performing calculations and setting values explicitly.
public class ResizableCanvas extends Pane {

    private final Canvas canvas = new Canvas();

    public ResizableCanvas() {
        getChildren().add(canvas);
        canvas.widthProperty().bind(this.prefWidthProperty());
        canvas.heightProperty().bind(this.prefHeightProperty());
    }

    public Canvas getCanvas() {
        return canvas;
    }
}
