package org.broad.igv.ui.javafx;

import javafx.scene.canvas.Canvas;
import javafx.scene.layout.Pane;

// Canvas doesn't resize on its own, so we need to wrap it in a container to handle this on its behalf
// http://fxexperience.com/2014/05/resizable-grid-using-canvas/
// https://stackoverflow.com/questions/27808210/java-fx-splitpane
// Warning: we need to watch out as there might be performance issues around this approach.
public class ResizableCanvas extends Pane {

  private final Canvas canvas = new Canvas();

  public ResizableCanvas() {
    getChildren().add(canvas);
  }

  @Override
  protected void layoutChildren() {
      final int top = (int)snappedTopInset();
      final int right = (int)snappedRightInset();
      final int bottom = (int)snappedBottomInset();
      final int left = (int)snappedLeftInset();
      final int w = (int)getWidth() - left - right;
      final int h = (int)getHeight() - top - bottom;
      canvas.setLayoutX(left);
      canvas.setLayoutY(top);
      if (w != canvas.getWidth() || h != canvas.getHeight()) {
          canvas.setWidth(w);
          canvas.setHeight(h);
      }
  }

  public Canvas getCanvas() {
    return canvas;
  }
}
