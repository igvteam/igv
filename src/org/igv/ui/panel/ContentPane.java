/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2018 Broad Institute
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
package org.igv.ui.panel;

import org.broad.igv.ui.panel.ReferenceFrame;
import org.igv.ui.JavaFXUIUtilities;
import org.igv.ui.ResizableCanvas;

/**
 * Superclass for Panes in the "Content" column of the UI (Cytobands, Rulers,
 * Tracks, etc.).
 *
 * @author eby
 */
public abstract class ContentPane extends ResizableCanvas {
    protected final ReferenceFrame frame;

    public ContentPane(ReferenceFrame frame) {
        this.frame = frame;
        JavaFXUIUtilities.bindWidthToProperty(this, frame.displayWidthProperty());
    }

    public ReferenceFrame getFrame() {
        return frame;
    }

    /**
     * Subclasses should call this as the final step of class construction; doing it
     * earlier may result in render() being called before subclass construction is
     * complete.
     *
     * Before exiting, this will make the initial render() call.
     */
    // Declared final to prevent overriding.
    protected final void completeInitialization() {
        // Re-render on change of width,height,zoom or chromosome.
        frame.chromosomeNameProperty().addListener((observable, oldValue, newValue) -> render());
        frame.zoomProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefWidthProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefHeightProperty().addListener((observable, oldValue, newValue) -> render());
        render();
    }

    abstract protected void render();
}
