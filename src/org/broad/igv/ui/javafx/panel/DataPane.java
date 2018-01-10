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

import javafx.event.EventHandler;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.input.MouseEvent;
import javafx.scene.paint.Color;

import org.apache.log4j.Logger;
import org.broad.igv.ui.javafx.JavaFXUIUtilities;
import org.broad.igv.ui.javafx.ResizableCanvas;
import org.broad.igv.ui.panel.ReferenceFrame;

// Intended as the rough equivalent of the DataPanel class of the Swing UI.  Work in progress.
public class DataPane extends ResizableCanvas {
    private static Logger log = Logger.getLogger(DataPane.class);
    
    private ReferenceFrame frame;
    private DataPaneContainer parent;

    public DataPane(ReferenceFrame frame, DataPaneContainer parent) {
        this.frame = frame;
        JavaFXUIUtilities.bindWidthToProperty(this, frame.displayWidthProperty());

        this.parent = parent;
        init();
        render();
        this.prefWidthProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefHeightProperty().addListener((observable, oldValue, newValue) -> render());
    }

    public ReferenceFrame getFrame() {
        return frame;
    }


    private void init() {

        this.setOnMouseClicked(new EventHandler<MouseEvent>() {
            @Override
            public void handle(MouseEvent event) {
                render();
            }
        });

    }

    public void render() {
        GraphicsContext gc = getCanvas().getGraphicsContext2D();

        gc.setFill(Color.BLUE);
        gc.fillRect(75,75,100,100);
    }




}
