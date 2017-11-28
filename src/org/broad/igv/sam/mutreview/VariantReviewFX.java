package org.broad.igv.sam.mutreview;

import javafx.embed.swing.JFXPanel;
import javafx.embed.swing.SwingFXUtils;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.control.Toggle;
import javafx.scene.control.ToggleGroup;
import javafx.scene.image.*;
import javafx.scene.layout.BorderPane;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;
import java.util.ResourceBundle;

public class VariantReviewFX {

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

        Toggle selectedButton = artifactGroup.getSelectedToggle();

        //metadata.score =
        //metatdata.scoreString =
        System.out.println();

    }


    public static void open(Frame parent, BufferedImage bufferedImage, VariantReviewMetadata metadata) throws IOException {

        SwingUtilities.invokeLater(() -> {
            try {

                JDialog frame = new JDialog(parent, "Preferences", true);
                JFXPanel fxPanel = new JFXPanel();
                FXMLLoader loader = new FXMLLoader(VariantReviewFX.class.getResource("VariantReview.fxml"));
                BorderPane pane = loader.load();

                VariantReviewFX controller = loader.getController();

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
                frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            } catch (IOException e) {
                e.printStackTrace();
            }

        });
    }

}
