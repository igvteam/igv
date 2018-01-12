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

package org.igv.ui.toolbar;

import javafx.scene.control.Slider;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

// As it stands this custom class may be unnecessary; all of this could be done
// with a stock instance modified by its owner.  May need to add to it later or
// further customize, though.
public class ZoomSlider extends Slider {

    private ReferenceFrame frame;

    public ZoomSlider() {
        this(null);
    }

    public ZoomSlider(ReferenceFrame frame) {
        this.frame = frame;
        this.setBlockIncrement(1.0);
        this.setShowTickMarks(true);
        this.setSnapToTicks(true);
        this.minProperty().set(0);
        this.maxProperty().bindBidirectional(getReferenceFrame().maxZoomProperty());
        this.valueProperty().bindBidirectional(getReferenceFrame().zoomProperty());
    }

    private ReferenceFrame getReferenceFrame() {
        if (frame == null) return FrameManager.getDefaultFrame();
        return frame;
    }
}
