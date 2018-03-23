package org.broad.igv.sam.mutreview;

import com.google.gson.Gson;
import javafx.embed.swing.JFXPanel;
import javafx.embed.swing.SwingFXUtils;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Toggle;
import javafx.scene.control.ToggleGroup;
import javafx.scene.image.*;
import javafx.scene.layout.BorderPane;
import org.apache.log4j.Logger;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;
import java.util.Map;
import java.util.ResourceBundle;

public class VariantReviewFX {

    private static Logger log = Logger.getLogger(VariantReviewFX.class);

    private JDialog dialog;

    private BufferedImage image;

    private VariantReviewMetadata metadata;

    @FXML // ResourceBundle that was given to the FXMLLoader
    private ResourceBundle resources;

    @FXML // URL location of the FXML file that was given to the FXMLLoader
    private URL location;

    @FXML // fx:id="artifactGroup"
    private ToggleGroup artifactGroup; // Value injected by FXMLLoader

    @FXML // fx:id="imageView"
    private ImageView imageView; // Value injected by FXMLLoader

    @FXML
        // This method is called by the FXMLLoader when initialization is complete
    void initialize() {
        assert artifactGroup != null : "fx:id=\"artifactGroup\" was not injected: check your FXML file 'VariantReview.fxml'.";
    }

    @FXML
    void submit() {

        RadioButton selectedButton = (RadioButton) artifactGroup.getSelectedToggle();

        if (selectedButton == null) {
            MessageUtils.showMessage("No call selected");
        } else {
            String selectedButtonText = selectedButton.getText().toLowerCase();
            int score = -1;
            if ("yes".equals(selectedButtonText)) {
                score = 0;
            } else if ("no".equals(selectedButtonText)) {
                score = 1;
            } else if ("unknown".equals(selectedButtonText)) {
                score = 2;
            }


            metadata.score = score;

            try {
                GoogleCloudStorageHelper.upload(image, metadata);

            } catch (IOException e) {

                MessageUtils.showErrorMessage("Error uploading data", e);
                log.error(e);
            }

            SwingUtilities.invokeLater(() -> {
                dialog.setVisible(false);
            });

        }
    }

    @FXML
    void cancel() {

        SwingUtilities.invokeLater(() -> dialog.setVisible(false));
    }


    public static void open(Frame parent, BufferedImage bufferedImage, VariantReviewMetadata metadata) throws IOException {

        SwingUtilities.invokeLater(() -> {
            try {

                JDialog frame = new JDialog(parent, "Preferences", true);
                JFXPanel fxPanel = new JFXPanel();
                FXMLLoader loader = new FXMLLoader(VariantReviewFX.class.getResource("VariantReview.fxml"));
                BorderPane pane = loader.load();

                VariantReviewFX controller = loader.getController();

                controller.dialog = frame;
                controller.image = bufferedImage;
                controller.metadata = metadata;

                javafx.scene.image.Image image = SwingFXUtils.toFXImage(bufferedImage, null);
                controller.imageView.setImage(image);


                fxPanel.setScene(new Scene(pane));

                frame.add(fxPanel);
                frame.pack();
                frame.setSize(950, 800);
                frame.setLocationRelativeTo(parent);
                frame.setVisible(true);
            } catch (IOException e) {
                e.printStackTrace();
            }

        });
    }

}
