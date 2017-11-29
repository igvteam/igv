package org.broad.igv.sam.mutreview;

import javafx.embed.swing.JFXPanel;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.image.ImageView;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;

/**
 * Created by jrobinso on 11/12/17.
 */
public class VariantReviewFX2 {


    public static void open(Frame parent) throws IOException {

        SwingUtilities.invokeLater(() -> {


            JDialog frame = new JDialog(parent, "Preferences", true);
            JFXPanel fxPanel = new JFXPanel();

            BorderPane pane = init();

            fxPanel.setScene(new Scene(pane));

            frame.add(fxPanel);
            frame.pack();
            frame.setSize(950, 800);
            frame.setLocationRelativeTo(parent);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);


        });
    }

    public static void main(String[] args) throws IOException {

        open(null);

    }


    static BorderPane init() {

        BorderPane borderPane = new BorderPane();
        borderPane.setPrefHeight(700.0);
        borderPane.setPrefWidth(950.0);
        borderPane.setPadding(new Insets(10, 10, 10, 10));

        HBox hBox = new HBox();
        borderPane.setTop(hBox);
        borderPane.setMargin(hBox, new Insets(10, 0, 0, 10));
        hBox.setPrefWidth(200.0);

        Label label = new Label("Variant is an artifact ?");
        label.setFont(new Font("System Bold", 24));
        hBox.setMargin(label, new Insets(0, 20, 0, 0));
        hBox.getChildren().add(label);

        VBox vbox = new VBox(10);
        hBox.getChildren().add(vbox);

        vbox.setPrefWidth(100);

        ToggleGroup group = new ToggleGroup();
        RadioButton yesButton = new RadioButton("Yes");
        yesButton.setToggleGroup(group);
        RadioButton noButton = new RadioButton("No");
        noButton.setToggleGroup(group);
        RadioButton unknownButton = new RadioButton("Unknown");
        unknownButton.setToggleGroup(group);

        vbox.getChildren().add(yesButton);
        vbox.getChildren().add(noButton);
        vbox.getChildren().add(unknownButton);


        ImageView imageView = new ImageView();
        imageView.setFitWidth(800);
        imageView.setFitHeight(500);
        borderPane.setCenter(imageView);

        ButtonBar buttonBar = new ButtonBar();
        buttonBar.setPadding(new Insets(10));
        borderPane.setBottom(buttonBar);

        buttonBar.setPrefWidth(200);
        buttonBar.setPrefHeight(40);

        Button cancelButton = new Button("Cancel");
        buttonBar.getButtons().add(cancelButton);

        Button sumbitButton = new Button("Submit");
        buttonBar.getButtons().add(sumbitButton);

        return borderPane;
    }
}
