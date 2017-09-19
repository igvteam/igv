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

import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.Session;
import org.broad.igv.ui.Main;
import org.broad.igv.ui.UIConstants;

import javafx.geometry.Rectangle2D;
import javafx.stage.Screen;
import javafx.stage.Stage;

// Intended as the rough equivalent of the IGV class of the Swing UI.  Work in progress.
// Here's a general list of remaining responsibilities from the origin class which may or may not end up 
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
public class IGVStageManager implements IGVEventObserver {
  
  private static IGVStageManager theInstance = null;
  
  // TODO: evaluate whether these can be final
  private Stage mainStage;
  private GenomeManager genomeManager;

  Session session;

  public static IGVStageManager createInstance(Stage stage) {
      if (theInstance != null) {
          throw new RuntimeException("Only a single instance is allowed.");
      }
      theInstance = new IGVStageManager(stage);
      return theInstance;
  }

  public static IGVStageManager getInstance() {
      if (theInstance == null) {
          throw new RuntimeException("IGVStageManager has not been initialized.  Must call createInstance(Stage) first");
      }
      return theInstance;
  }

  public Stage getMainStage() {
    return mainStage;
  }

  private IGVStageManager(Stage stage) {
    mainStage = stage;
    
    // TODO: port this to JavaFX (AWT refs)
    final IGVPreferences preferences = PreferencesManager.getPreferences();

    // TODO: port this to JavaFX (Swing refs)
    genomeManager = GenomeManager.getInstance();
    
    // TODO: add JavaFX equiv. of WindowListener for handling loss of focus
    
    // TODO: review this for porting to JavaFX 
    session = new Session(null);

    // TODO: necessary to create cursors in JavaFX?  Maybe for easy reference from elsewhere.  Do it later.
    
    mainStage.setTitle(UIConstants.APPLICATION_NAME);
    
    // TODO: Need equivalent of GlassPane for the Scene
    // Note: we can use an EventFilter on the Scene instead of doing it as a GlassPane.  Suggested here:
    // https://stackoverflow.com/questions/18758803/javafx-best-practice-to-implement-glasspane-like-mouse-touch-event-receiver-whi
    IGVContentManager igvContentManager = new IGVContentManager(this);
    
    // Certain components MUST be visible, so we set minimum size
    // TODO: evaluate if this is the same for JavaFX
    mainStage.setMinHeight(300);
    mainStage.setMinWidth(300);
    
    // Set the application's previous location and size; 
    Rectangle2D screenBounds = Screen.getPrimary().getVisualBounds();
    Rectangle2D applicationBounds = preferences.getApplicationFrameBounds_javafx();

    if (applicationBounds == null || applicationBounds.getMaxX() > screenBounds.getWidth() ||
            applicationBounds.getMaxY() > screenBounds.getHeight() ||
            applicationBounds.getWidth() == 0 || applicationBounds.getHeight() == 0) {
        int width = Math.min(1150, (int) screenBounds.getWidth());
        int height = Math.min(800, (int) screenBounds.getHeight());
        applicationBounds = new Rectangle2D(0, 0, width, height);
    }
    mainStage.setX(applicationBounds.getMaxX());
    mainStage.setY(applicationBounds.getMaxY());
    mainStage.setWidth(applicationBounds.getWidth());
    mainStage.setHeight(applicationBounds.getHeight());
  }
  
  public void startUp(Main.IGVArgs igvArgs) {
    // TODO: implement this in terms of JavaFX concurrency (Task? Worker?).
    // Just swallowing this call for now.
  }
  
  @Override
  public void receiveEvent(Object event) {
    // TODO Auto-generated method stub
    
  }
}
