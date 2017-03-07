package org.broad.igv.prefs;

import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.embed.swing.JFXPanel;
import javafx.geometry.HPos;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TextField;
import javafx.scene.layout.*;
import javafx.stage.Stage;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.beans.EventHandler;
import java.io.*;
import java.util.*;
import java.util.List;

public class PreferenceEditorFX {

    public static void main(String[] args) throws IOException {

        open(null);

    }

    public static void open(Frame parent) throws IOException {
        List<PreferencesManager.PreferenceGroup> preferenceGroups = PreferencesManager.loadPreferenceList();

        SwingUtilities.invokeLater(() -> {
            JDialog frame = new JDialog(parent, "Preferences", true);

            final JFXPanel fxPanel = new JFXPanel();
            Platform.runLater(() -> initFX(frame, fxPanel, preferenceGroups));
            frame.add(fxPanel);
            frame.pack();
            frame.setSize(800, 600);
            frame.setLocationRelativeTo(parent);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        });
    }

    private static void initFX(final JDialog parent, JFXPanel fxPanel, List<PreferencesManager.PreferenceGroup> preferenceGroups) {

        // final Map<String, String> updatedPrefs = new HashMap<>();

        TabPane tabPane = new TabPane();
        BorderPane borderPane = new BorderPane();
        borderPane.setCenter(tabPane);
        Scene scene = new Scene(borderPane);

        Map<String, Map<String, String>> updatedPreferencesMap = new HashMap<>();
        for (PreferencesManager.PreferenceGroup entry : preferenceGroups) {

            if(entry.tabLabel.equals("Hidden")) continue;

            final IGVPreferences preferences = PreferencesManager.getPreferences(entry.category);

            final Map<String, String> updatedPrefs =
                    updatedPreferencesMap.containsKey(entry.category) ? updatedPreferencesMap.get(entry.category) :
                            new HashMap<>();
            updatedPreferencesMap.put(entry.category, updatedPrefs);


            final String tabLabel = entry.tabLabel;
            Tab tab = new Tab(tabLabel);
            tab.setClosable(false);

            tabPane.getTabs().add(tab);

            ScrollPane scrollPane = new ScrollPane();
            tab.setContent(scrollPane);

            VBox vBox = new VBox();
            vBox.setFillWidth(true);
            //   vBox.prefWidthProperty().bind(pane.widthProperty());
            scrollPane.setContent(vBox);

            GridPane gridPane = new GridPane();
            gridPane.setHgap(5);
            gridPane.setVgap(5);
            vBox.getChildren().add(gridPane);

            String currentGroup = null;

            int row = 1;
            for (PreferencesManager.Preference pref : entry.preferences) {

                try {
                    Tooltip tooltip = pref.getComment() == null ? null : new Tooltip(pref.getComment());

                    if (pref.getKey().equals(PreferencesManager.SEPARATOR_KEY)) {
                        Separator sep = new Separator();
                        GridPane.setColumnSpan(sep, 4);
                        gridPane.add(sep, 1, row);
                        row++;
                        continue;
                    }

                    if (pref.getKey().equals(PreferencesManager.INFO_KEY)) {
                        row++;
                        Label label = new Label(pref.getLabel());
                        label.setStyle("-fx-font-size:16");
                        label.setStyle("-fx-font-weight: bold");
                        GridPane.setColumnSpan(label, 4);
                        gridPane.add(label, 1, row);
                        row += 2;
                        continue;
                    }

                    if (pref.group != null && !pref.group.equals(currentGroup)) {
                        // Start a new group
                        row = 0;

                        currentGroup = pref.group;
                        gridPane = new GridPane();
                        gridPane.setHgap(5);
                        gridPane.setVgap(5);

                        vBox.getChildren().add(gridPane);
                        TitledPane tp = new TitledPane(currentGroup, gridPane);
                        tp.setCollapsible(false);

                        vBox.getChildren().add(tp);
                        vBox.setMargin(tp, new Insets(10, 0, 0, 0));

                    }

                    if (pref.getType().equals("boolean")) {
                        CheckBox cb = new CheckBox(pref.getLabel());
                        cb.setSelected(preferences.getAsBoolean(pref.getKey()));
                        cb.setOnAction(event -> {
                            updatedPrefs.put(pref.getKey(), Boolean.toString(cb.isSelected()));
                            System.out.println("Set " + pref.getLabel() + ": " + cb.isSelected());
                        });
                        GridPane.setColumnSpan(cb, 2);
                        gridPane.add(cb, 1, row);

                        if (tooltip != null) {
                            cb.setTooltip(tooltip);
                        }

                    } else if (pref.getType().startsWith("select")) {
                        Label label = new Label(pref.getLabel());
                        String[] selections = Globals.whitespacePattern.split(pref.getType())[1].split("\\|");
                        final ComboBox comboBox = new ComboBox(FXCollections.observableArrayList(Arrays.asList(selections)));
                        comboBox.valueProperty().setValue(pref.getDefaultValue());
                        comboBox.valueProperty().addListener((ov, t, t1) -> {
                            System.out.println("Set " + pref.getLabel() + " " + comboBox.valueProperty().toString());
                        });
                        gridPane.add(label, 1, row);
                        GridPane.setColumnSpan(comboBox, 3);
                        gridPane.add(comboBox, 2, row);

                        if (tooltip != null) {
                            label.setTooltip(tooltip);
                            comboBox.setTooltip(tooltip);
                        }

                    } else {
                        Label label = new Label(pref.getLabel());
                        TextField field = new TextField(preferences.get(pref.getKey()));
                        field.setPrefWidth(500);
                        field.setOnAction(event -> {
                            final String text = field.getText();
                            if (validate(text, pref.getType())) {
                                updatedPrefs.put(pref.getKey(), text);
                            } else {
                                field.setText(preferences.get(pref.getKey()));
                            }
                        });
                        field.focusedProperty().addListener((observable, oldValue, newValue) -> {
                            if (newValue == false) {
                                final String text = field.getText();
                                if (validate(text, pref.getType())) {
                                    updatedPrefs.put(pref.getKey(), text);
                                } else {
                                    field.setText(preferences.get(pref.getKey()));
                                }
                            }
                        });

                        gridPane.add(label, 1, row);
                        gridPane.add(field, 2, row);


                        if (tooltip != null) {
                            label.setTooltip(tooltip);
                            field.setTooltip(tooltip);
                        }
                    }


                    row++;
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            if (tabLabel.equalsIgnoreCase("Advanced")) {
                // Add IGV directory management at the end.  This is a special case
                String currentDirectory = DirectoryManager.getIgvDirectory().getAbsolutePath();
                final Label currentDirectoryLabel = new Label("IGV Directory: " + currentDirectory);
                final Button moveButton = new Button("Move...");
                row++;
                gridPane.add(currentDirectoryLabel, 1, row);
                GridPane.setHalignment(moveButton, HPos.RIGHT);
                gridPane.add(moveButton, 2, row);

                moveButton.setOnAction(event -> {
                    // Do this on the Swing thread until we port to javafx file dialog
                    UIUtilities.invokeOnEventThread(() -> {
                        final File igvDirectory = DirectoryManager.getIgvDirectory();
                        final File newDirectory = FileDialogUtils.chooseDirectory("Select IGV directory", DirectoryManager.getUserDirectory());
                        if (newDirectory != null && !newDirectory.equals(igvDirectory)) {
                            DirectoryManager.moveIGVDirectory(newDirectory);
                            Platform.runLater(() -> currentDirectoryLabel.setText(newDirectory.getAbsolutePath()));
                        }
                    });
                });

            }
        }

        HBox hbox = new HBox();
        hbox.setAlignment(Pos.CENTER_RIGHT);
        hbox.setPadding(new Insets(15, 12, 15, 12));
        hbox.setSpacing(5);
        hbox.setStyle("-fx-background-color: #336699;");

        Button cancelBUtton = new Button("Cancel");
        cancelBUtton.setPrefSize(100, 20);
        cancelBUtton.setOnAction((event) -> {
            SwingUtilities.invokeLater(() -> parent.setVisible(false));
        });

        Button saveButton = new Button("Save");
        saveButton.setPrefSize(100, 20);
        saveButton.setDefaultButton(true);
        saveButton.setOnAction((event) -> {
            PreferencesManager.updateAll(updatedPreferencesMap);
            SwingUtilities.invokeLater(() -> parent.setVisible(false));
            if(IGV.hasInstance()) {
                IGV.getInstance().doRefresh();
            }
        });
        hbox.getChildren().addAll(cancelBUtton, saveButton);

        borderPane.setBottom(hbox);

        fxPanel.setScene(scene);
    }


    private static boolean validate(String text, String type) {

        if (type.equals("integer")) {
            try {
                Integer.parseInt(text);
            } catch (NumberFormatException e) {
                return false;
            }
        } else if (type.equals("float")) {
            try {
                Double.parseDouble(text);
            } catch (NumberFormatException e) {
                return false;
            }
        }

        return true;
    }


}
