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
package org.broad.igv.ui.javafx.toolbar;

import javafx.application.Platform;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.scene.control.ComboBox;
import javafx.scene.control.ListCell;
import javafx.scene.control.ListView;
import javafx.scene.text.Text;
import javafx.util.Callback;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.javafx.FontMetrics;
import org.broad.igv.ui.panel.FrameManager;

import java.util.List;

/**
 * @author jrobinso on 7/6/17.
 * @author eby rewritten for JavaFX
 */
public class ChromosomeComboBox extends ComboBox<String> {
    private static Logger log = Logger.getLogger(ChromosomeComboBox.class);
    private final static double DEFAULT_CHROMOSOME_DROPDOWN_WIDTH = 120;

    private DoubleProperty controlWidthProperty = new SimpleDoubleProperty(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH);

    public ChromosomeComboBox(Genome genome) {
        prefWidthProperty().bind(controlWidthProperty);
        minWidthProperty().bind(controlWidthProperty);
        maxWidthProperty().bind(controlWidthProperty);

        valueProperty().bindBidirectional(FrameManager.getDefaultFrame().chromosomeNameProperty());
        updateChromosFromGenome(genome);
    }

    // Should be able to replace this with binding to an ObjectProperty<Genome>
    public void updateChromosFromGenome(Genome genome) {

        if (genome == null) return;

        Platform.runLater(() -> {
            // Reset the whole control to default width
            controlWidthProperty.set(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH);
            double w = DEFAULT_CHROMOSOME_DROPDOWN_WIDTH;

            List<String> allChromosomeNames = genome.getAllChromosomeNames();

            setCellFactory(new Callback<ListView<String>, ListCell<String>>() {

                @Override
                public ListCell<String> call(ListView<String> param) {
                    return new ListCell<String>() {

                        @Override
                        protected void updateItem(String item, boolean empty) {
                            super.updateItem(item, empty);
                            if (!empty && item != null) {
                                setText(item);
                            }
                        }
                    };
                }
            });

            getItems().clear();
            getItems().addAll(allChromosomeNames);
            if (!allChromosomeNames.isEmpty()) {
                String homeChr = genome.getHomeChromosome();
                if (homeChr.equals(Globals.CHR_ALL)) {
                    getItems().add(0, Globals.CHR_ALL);
                }
            }
            getSelectionModel().select(genome.getHomeChromosome());

            Text sizer = new Text("");
            sizer.setFont(getEditor().getFont());
            for (String item : allChromosomeNames) {
                double width = FontMetrics.getTextWidthInFont(item, sizer) + 50;
                if (width > w) {
                    w = width;
                }
            }
            if (w > DEFAULT_CHROMOSOME_DROPDOWN_WIDTH) {
                controlWidthProperty.set(w);
            }
        });
    }
}
