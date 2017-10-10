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

import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.geometry.Orientation;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.SplitPane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.StrokeLineCap;
import javafx.stage.Screen;
import javafx.stage.Stage;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.javafx.panel.MainContentPane;

// Builds the Stage content.  This was originally intended as the rough equivalent of the IGV class of the Swing UI,
// but for now is just populating the Stage.  In the long run the other IGV-equivalent responsibilities may land here,
// but that will become more clear as we progress.
// Work in progress. Here's a general list of remaining responsibilities from the origin class which may or may not end up 
// being re-implemented here (they will go somewhere, maybe elsewhere):
// TODO: event handling for toolbar items
// TODO: event handling for GlassPane (what ever that does).  See notes below.
// TODO: event handling for various menu items
// TODO: cursor changes; wait, zoom, hand, DnD
// TODO: Drag and Drop handling (for ROI, I think)
// TODO: manage overlaid tracks
// TODO: handle save Prefs on exit
// TODO: populate controls (e.g. Genome & Chromosome drop-downs)
// TODO: manage Tracks
// TODO: launch initial Task based on the cmdLine params
public class IGVStageBuilder {

    public static void buildStage(Stage stage) {
        // TODO: port this to JavaFX (AWT refs)
        final IGVPreferences preferences = PreferencesManager.getPreferences();

        // TODO: add JavaFX equiv. of WindowListener for handling loss of focus

        // TODO: necessary to create cursors in JavaFX? Maybe for easy reference from
        // elsewhere. Do it later.

        stage.setTitle(UIConstants.APPLICATION_NAME);

        // Certain components MUST be visible, so we set minimum size
        // TODO: evaluate if this is the same for JavaFX
        stage.setMinHeight(300);
        stage.setMinWidth(300);

        // Set the application's previous location and size;
        Rectangle2D screenBounds = Screen.getPrimary().getVisualBounds();
        Rectangle2D applicationBounds = preferences.getApplicationFrameBounds_javafx();

        if (applicationBounds == null || applicationBounds.getMaxX() > screenBounds.getWidth()
                || applicationBounds.getMaxY() > screenBounds.getHeight() || applicationBounds.getWidth() == 0
                || applicationBounds.getHeight() == 0) {
            int width = Math.min(1150, (int) screenBounds.getWidth());
            int height = Math.min(800, (int) screenBounds.getHeight());
            applicationBounds = new Rectangle2D(0, 0, width, height);
        }
        stage.setX(applicationBounds.getMaxX());
        stage.setY(applicationBounds.getMaxY());
        stage.setWidth(applicationBounds.getWidth());
        stage.setHeight(applicationBounds.getHeight());

        // TODO: Need equivalent of GlassPane for the Scene
        // Note: we can use an EventFilter on the Scene instead of doing it as a
        // GlassPane. Suggested here:
        // https://stackoverflow.com/questions/18758803/javafx-best-practice-to-implement-glasspane-like-mouse-touch-event-receiver-whi
        Scene contentScene = buildContent(stage);
        stage.setScene(contentScene);
    }

    private static Scene buildContent(Stage stage) {

        VBox contentContainer = new VBox();
        Scene contentScene = new Scene(contentContainer, stage.getHeight(), stage.getWidth(), Color.WHITE);
        contentScene.getStylesheets().add(IGVStageBuilder.class.getResource("igv.css").toExternalForm());

        IGVMenuBarManager igvMenuBarBuilder = new IGVMenuBarManager(stage);
        contentContainer.getChildren().add(igvMenuBarBuilder.getMenuBar());

        IGVToolBarManager igvToolBar = new IGVToolBarManager();
        contentContainer.getChildren().add(igvToolBar.getToolBar());

        MainContentPane mainContentPane = new MainContentPane();
        contentContainer.getChildren().add(mainContentPane);
        VBox.setVgrow(mainContentPane, Priority.ALWAYS);
        mainContentPane.prefWidthProperty().bind(contentContainer.widthProperty());
        // Replace with actual content grid
        //drawPlaceholderContent(contentContainer);

        stage.showingProperty().addListener(new ChangeListener<Boolean>() {

            @Override
            public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                if (newValue) {
                    mainContentPane.addInitialTrackRows();

                    observable.removeListener(this);
                }
            }
        });

        return contentScene;
    }

    private static void drawPlaceholderContent(VBox mainPanel) {
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
