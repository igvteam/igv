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
import javafx.scene.control.*;
import javafx.scene.layout.HBox;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.igv.event.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.javafx.toolbar.ChromosomeComboBox;
import org.broad.igv.ui.javafx.toolbar.ZoomSlider;
import org.broad.igv.ui.panel.FrameManager;

// Intended as the rough equivalent of the IGVCommandBar class of the Swing UI.  Work in progress.
// Will add event handlers (or at least stubs) for all of the included controls.
public class IGVToolBarManager implements IGVEventObserver {
    private static Logger log = Logger.getLogger(IGVToolBarManager.class);
    
    private ToolBar toolBar;
    private ComboBox<String> genomeSelector;
    private ChromosomeComboBox chromosomeSelector;
    private ZoomSlider zoomSlider;

    public IGVToolBarManager() {

        // TODO: populate ToolBar controls with actual content.
        // TODO: add event handlers to all ToolBar components, including enable/disable.
        // Probably need to migrate most controls to instance variables.
        String defaultGenome = "Human hg19";
        ObservableList<String> genomes = FXCollections.observableArrayList(defaultGenome, "chr1.fasta");
        genomeSelector = new ComboBox<String>(genomes);

        chromosomeSelector = new ChromosomeComboBox(GenomeManager.getInstance().getCurrentGenome());
        
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

        zoomSlider = new ZoomSlider();
        zoomSlider.setPrefWidth(200.0);
        toolBar = new ToolBar(genomeSelector, chromosomeSelector, jumpToPane,
                homeButton, leftArrowButton, rightArrowButton, refreshScreenButton, regionToolButton, resizeToWindowButton,
                infoSelectButton, rulerButton, zoomSlider);

        log.info("About to register with eventBus");
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(GenomeResetEvent.class, this);
        log.info("Done registering with eventBus");
    }

    public ToolBar getToolBar() {
        return toolBar;
    }

    @Override
    public void receiveEvent(Object e) {
        if (e instanceof ViewChange) {
            ViewChange event = (ViewChange) e;
            log.info("View Change event: " + event.type);
            if (event.type == ViewChange.Type.ChromosomeChange || event.type == ViewChange.Type.LocusChange) {
                String chrName = FrameManager.getDefaultFrame().getChrName();
                log.info("change chr to: " + chrName);
                // TODO: enable ROI toggle & zoom control
                if (StringUtils.isNotBlank(chrName) && !chrName.equals(chromosomeSelector.getSelectionModel().getSelectedItem())) {
                    log.info("change chrSelector to: " + chrName);
                    chromosomeSelector.getSelectionModel().select(chrName);
                }

                updateCurrentCoordinates();
            }
        } else if (e instanceof GenomeChangeEvent) {
            GenomeChangeEvent event = (GenomeChangeEvent) e;
            Genome genome = event.genome;
            refreshGenomeListComboBox();
            chromosomeSelector.updateChromosFromGenome(genome);
        } else if (e instanceof GenomeResetEvent) {
            refreshGenomeListComboBox();
        } else {
            log.info("Unknown event class: " + e.getClass());
        }
    }

    public void refreshGenomeListComboBox() {
        // Not yet implemented
    }

    public void updateCurrentCoordinates() {
        // Not yet implemented
    }
}
