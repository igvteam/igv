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

import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.control.ToolBar;
import javafx.scene.layout.HBox;

// Intended as the rough equivalent of the IGVCommandBar class of the Swing UI.  Work in progress.
public class IGVToolBarManager {

  private ToolBar toolBar;
  
  public IGVToolBarManager() {

    // TODO: populate ToolBar controls with actual content.
    // TODO: add event handlers to all ToolBar components, including enable/disable.
    // Do these need to be held as instance vars for better access from event handlers?  Not sure of the 
    // best structure here yet.
    String defaultGenome = "Human hg19";
    ObservableList<String> genomes = FXCollections.observableArrayList(defaultGenome, "chr1.fasta");
    ComboBox<String> genomeSelector = new ComboBox<String>(genomes);
    genomeSelector.setValue(defaultGenome);

    ObservableList<String> chromosomes = FXCollections.observableArrayList("chr1, chr2, chr3");
    ComboBox<String> chromosomeSelector = new ComboBox<String>(chromosomes);
    
    TextField jumpToTextField = new TextField();
    Label jumpToLabel = new Label("Go");
    jumpToLabel.setLabelFor(jumpToTextField);
    HBox jumpToPane = new HBox(jumpToTextField, jumpToLabel);
    jumpToPane.setAlignment(Pos.CENTER);

    Button homeButton = new Button("");
    homeButton.setId("homeButton");
    Button leftArrowButton = new Button("");
    leftArrowButton.setId("leftArrowButton");
    Button rightArrowButton = new Button("");
    rightArrowButton.setId("rightArrowButton");
    Button refreshScreenButton = new Button("");
    refreshScreenButton.setId("refreshScreenButton");
    Button regionToolButton = new Button("");
    regionToolButton.setId("regionToolButton");
    Button resizeToWindowButton = new Button("");
    resizeToWindowButton.setId("resizeToWindowButton");
    Button infoSelectButton = new Button("");
    infoSelectButton.setId("infoSelectButton");
    Button rulerButton = new Button("");
    rulerButton.setId("rulerButton");
    
    Slider zoomLevelSlider = new Slider();
    zoomLevelSlider.setShowTickMarks(true);
    zoomLevelSlider.setSnapToTicks(true);

    toolBar = new ToolBar(genomeSelector, chromosomeSelector, jumpToPane,
        homeButton, leftArrowButton, rightArrowButton, refreshScreenButton, regionToolButton, resizeToWindowButton,
        infoSelectButton, rulerButton, zoomLevelSlider);
  }

  public ToolBar getToolBar() {
    return toolBar;
  }
}
