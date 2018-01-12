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

import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.ToolBar;
import javafx.scene.layout.HBox;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.igv.ui.toolbar.ChromosomeComboBox;
import org.igv.ui.toolbar.SearchTextField;
import org.igv.ui.toolbar.ZoomSlider;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

// Intended as the rough equivalent of the IGVCommandBar class of the Swing UI.  Work in progress.
// Will add event handlers (or at least stubs) for all of the included controls.
public class IGVToolBarManager implements IGVEventObserver {
    private static Logger log = Logger.getLogger(IGVToolBarManager.class);
    
    private ToolBar toolBar;
    private ComboBox<String> genomeSelector;
    private ChromosomeComboBox chromosomeSelector;
    private ZoomSlider zoomSlider;
    private SearchTextField searchTextField = new SearchTextField();
    private Button goButton = new Button("Go");

    public IGVToolBarManager() {

        // TODO: populate ToolBar controls with actual content.
        // TODO: add event handlers to all ToolBar components, including enable/disable.
        // Probably need to migrate most controls to instance variables.
        String defaultGenome = "Human hg19";
        ObservableList<String> genomes = FXCollections.observableArrayList(defaultGenome, "chr1.fasta");
        genomeSelector = new ComboBox<String>(genomes);

        chromosomeSelector = new ChromosomeComboBox(GenomeManager.getInstance().getCurrentGenome());

        goButton.setOnAction((event) -> searchByLocus(searchTextField.getText()));
        goButton.setId("goButton");
        HBox jumpToPane = new HBox(3, searchTextField, goButton);
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

        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(GenomeResetEvent.class, this);
    }

    public ToolBar getToolBar() {
        return toolBar;
    }

    @Override
    public void receiveEvent(Object e) {
        if (e instanceof ViewChange) {
            ViewChange event = (ViewChange) e;
            if (event.type == ViewChange.Type.ChromosomeChange || event.type == ViewChange.Type.LocusChange) {
                String chrName = FrameManager.getDefaultFrame().getChrName();
                // TODO: enable ROI toggle & zoom control
                if (StringUtils.isNotBlank(chrName) && !chrName.equals(chromosomeSelector.getSelectionModel().getSelectedItem())) {
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
        String p = "";

        ReferenceFrame defaultFrame = FrameManager.getDefaultFrame();
        final String chrName = defaultFrame.getChrName();
        if (!Globals.CHR_ALL.equals(chrName) && !FrameManager.isGeneListMode()) {
            p = defaultFrame.getFormattedLocusString();
        }
        final String position = p;
        // Not dealing with history or IGV singleton yet
        //final History history = IGV.getInstance().getSession().getHistory();
        Platform.runLater(() -> {
            searchTextField.setText(position);
            //forwardButton.setEnabled(history.canGoForward());
            //backButton.setEnabled(history.canGoBack());
            //roiToggleButton.setEnabled(!Globals.CHR_ALL.equals(chrName));
            //zoomControl.setEnabled(!Globals.CHR_ALL.equals(chrName));
        });
    }

    public void searchByLocus(final String searchText) {

        if ((searchText != null) && (searchText.length() > 0)) {
            String homeChr = FrameManager.getDefaultFrame().getChrName();
            if (searchText.equalsIgnoreCase("home") || searchText.equalsIgnoreCase(homeChr)) {
                // Not yet implemented
                //homeButtonActionPerformed(null);
            } else {
                searchTextField.setText(searchText);
                searchTextField.searchByLocus(searchText);
            }
        }
    }
}
