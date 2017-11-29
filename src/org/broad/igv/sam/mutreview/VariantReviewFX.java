package org.broad.igv.sam.mutreview;

import javafx.embed.swing.JFXPanel;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.control.Toggle;
import javafx.scene.control.ToggleGroup;
import javafx.scene.layout.BorderPane;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.net.URL;
import java.util.ResourceBundle;

public class VariantReviewFX {

    @FXML // ResourceBundle that was given to the FXMLLoader
    private ResourceBundle resources;

    @FXML // URL location of the FXML file that was given to the FXMLLoader
    private URL location;

    @FXML // fx:id="artifactGroup"
    private ToggleGroup artifactGroup; // Value injected by FXMLLoader

    @FXML // This method is called by the FXMLLoader when initialization is complete
    void initialize() {
        assert artifactGroup != null : "fx:id=\"artifactGroup\" was not injected: check your FXML file 'VariantReview.fxml'.";
    }

    @FXML
    void submit() {
        Toggle selectedButton = artifactGroup.getSelectedToggle();

    }


    public static void open(Frame parent) throws IOException {

        SwingUtilities.invokeLater(() -> {
            try {

                JDialog frame = new JDialog(parent, "Preferences", true);
                JFXPanel fxPanel = new JFXPanel();
                BorderPane pane = FXMLLoader.load(VariantReviewFX.class.getResource("VariantReview.fxml"));

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

    public static void main(String[] args) throws IOException {

        open(null);

    }
}
