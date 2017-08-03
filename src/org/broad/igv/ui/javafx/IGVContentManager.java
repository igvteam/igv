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

import javafx.geometry.Orientation;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.SplitPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.StrokeLineCap;

// Intended as the rough equivalent of the IGVContentPane class of the Swing UI.  Work in progress.
public class IGVContentManager {
  
  // Since this is a Singleton, is it necessary to hold a reference?
  //private IGVStageManager igvStageManager;
  
  private Scene contentScene;
  
  public IGVContentManager(IGVStageManager igvStageManager) {
    //this.igvStageManager = igvStageManager;

    VBox contentContainer = new VBox();
    contentScene = new Scene(contentContainer, 300, 300, Color.WHITE);
    contentScene.getStylesheets().add(IGVContentManager.class.getResource("igv.css").toExternalForm());
    
    IGVMenuBarManager igvMenuBarManager = IGVMenuBarManager.createInstance(igvStageManager);
    contentContainer.getChildren().add(igvMenuBarManager.getMenuBar());
    
    IGVToolBarManager igvToolBarManager = new IGVToolBarManager();
    contentContainer.getChildren().add(igvToolBarManager.getToolBar());
    
    drawPlaceholderContent(contentContainer);
    
    igvStageManager.getMainStage().setScene(contentScene);
  }
  
  public Scene getContentScene() {
    return contentScene;
  }

  private void drawPlaceholderContent(VBox mainPanel) {
    // TODO: move the canvases to dedicated components dealing with real content
    ResizableCanvas resizableCanvas1 = new ResizableCanvas();
    Canvas canvas = resizableCanvas1.getCanvas();
    canvas.setHeight(500);
    canvas.setWidth(1000);
    
    ResizableCanvas resizableCanvas2 = new ResizableCanvas();
    Canvas canvas2 = resizableCanvas2.getCanvas();
    canvas2.setHeight(200);
    canvas2.setWidth(1000);
    
    SplitPane trackPane = new SplitPane(resizableCanvas1, resizableCanvas2);
    trackPane.setOrientation(Orientation.VERTICAL);
    trackPane.setDividerPositions(0.9);
    mainPanel.getChildren().add(trackPane);
    
    // Drawing placeholder content for now, nothing real
    GraphicsContext graphicsContext = canvas.getGraphicsContext2D();
    
    graphicsContext.setStroke(Color.RED);
    graphicsContext.setLineWidth(1);
    graphicsContext.setLineCap(StrokeLineCap.BUTT);
    graphicsContext.setLineDashes(10d, 5d, 15d, 20d);
    graphicsContext.setLineDashOffset(0);
    graphicsContext.strokeLine(10, 10, 200, 10);
    
    graphicsContext.setStroke(Color.GREEN);
    graphicsContext.setLineWidth(1);
    graphicsContext.setLineCap(StrokeLineCap.ROUND);
    graphicsContext.strokeLine(10, 30, 200, 30);
    
    graphicsContext.setStroke(Color.BLUE);
    graphicsContext.setLineWidth(1);
    graphicsContext.strokeLine(10, 50, 200, 50);

    GraphicsContext graphicsContext2 = canvas2.getGraphicsContext2D();
    
    graphicsContext2.setStroke(Color.PURPLE);
    graphicsContext2.setLineWidth(1);
    graphicsContext2.setLineCap(StrokeLineCap.BUTT);
    graphicsContext2.strokeLine(10, 10, 200, 10);
  }

}
