package org.broad.igv.ui.javafx.toolbar;

import javafx.application.Platform;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.event.ActionEvent;
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

    public ChromosomeComboBox() {
        setOnAction(evt -> chromosomeComboBoxActionPerformed(evt));

        setPrefHeight(16);
        setMinHeight(27);
        setMaxHeight(35);

        prefWidthProperty().bind(controlWidthProperty);
        minWidthProperty().bind(controlWidthProperty);
        maxWidthProperty().bind(controlWidthProperty);
    }

    private void chromosomeComboBoxActionPerformed(ActionEvent evt) {
        final String chrName = this.getSelectionModel().getSelectedItem().toString();
        if (chrName != null & !chrName.equals(FrameManager.getDefaultFrame().getChrName())) {
            FrameManager.getDefaultFrame().changeChromosome(chrName, true);
        }
    }

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

    class AutoWidthResizeListCell extends ListCell<String> {
        public AutoWidthResizeListCell(List<String> items, DoubleProperty controlWidthProperty, double defaultWidth) {
            // Reset the whole control to default width
            controlWidthProperty.set(defaultWidth);
            double w = defaultWidth;

            // Bind this ListCell's width to that of the parent control
            prefWidthProperty().bind(controlWidthProperty);
            maxWidthProperty().bind(controlWidthProperty);
            minWidthProperty().bind(controlWidthProperty);

            Text sizer = new Text("");
            sizer.setFont(getFont());
            for (String item : items) {
                double width = FontMetrics.getTextWidthInFont(item, sizer) + 50;
                if (width > w) {
                    w = width;
                    controlWidthProperty.set(w);
                }
            }
        }

        @Override
        protected void updateItem(String item, boolean empty) {
            super.updateItem(item, empty);
            if (!empty && item != null) {
                setText(item);
            }
        }


    }
}
