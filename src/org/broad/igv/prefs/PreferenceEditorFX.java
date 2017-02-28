package org.broad.igv.prefs;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.embed.swing.JFXPanel;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import oracle.jdbc.proxy.annotation.Pre;
import org.apache.batik.util.PreferenceManager;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class PreferenceEditorFX {

    public static final String SEPARATOR_KEY = "---";

    public static void main(String[] args) throws IOException {

        Map<String, List<Preference>> prefMap = loadPreferenceList();

        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Preferences");
            final JFXPanel fxPanel = new JFXPanel();
            frame.add(fxPanel);
            frame.setSize(800, 600);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

            Platform.runLater(() -> initFX(fxPanel, prefMap));
        });
    }

    private static void initFX(JFXPanel fxPanel, Map<String, List<Preference>> prefMap) {

        final IGVPreferences preferences = PreferencesManager.getPreferences();
        final Map<String, String> updatedPrefs = new HashMap<>();

        TabPane pane = new TabPane();
        Scene scene = new Scene(pane);

        for (Map.Entry<String, List<Preference>> entry : prefMap.entrySet()) {

            final String tabLabel = entry.getKey();
            Tab tab = new Tab(tabLabel);
            tab.setClosable(false);

            pane.getTabs().add(tab);

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
            for (Preference pref : entry.getValue()) {

                try {
                    if (pref.getKey().equals(SEPARATOR_KEY)) {
                        Separator sep = new Separator();
                        GridPane.setColumnSpan(sep, 4);
                        gridPane.add(sep, 1, row);
                        row++;
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

                    } else {
                        Label label = new Label(pref.getLabel());
                        TextField field = new TextField(preferences.get(pref.getKey()));
                        field.setOnAction(event -> {
                            final String text = field.getText();
                            if (validate(text, pref.getType())) {
                                System.out.println("Set " + pref.getLabel() + ": " + text);
                                updatedPrefs.put(pref.getKey(), text);
                            } else {
                                field.setText(preferences.get(pref.getKey()));
                            }

                        });
                        gridPane.add(label, 1, row);
                        gridPane.add(field, 2, row);
                    }

                    if (pref.getComment() != null) {
                        Text label = new Text(pref.getComment());
                        label.setStyle("-fx-font-style: italic;");  // Currently doesn't work (JavaFX bug)
                        gridPane.add(label, 3, row);
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
                GridPane.setColumnSpan(currentDirectoryLabel, 3);
                gridPane.add(currentDirectoryLabel, 1, row);
                gridPane.add(moveButton, 3, row);

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


    private static class Preference {

        int tokenCount = 0;
        String group;
        String key;
        String[] tokens;


        Preference(String key, String group) {
            this.key = key;
            tokens = new String[4];
            this.group = group;
        }

        void pushToken(String token) {
            try {
                tokens[tokenCount++] = token;
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        String getKey() {
            return key;
        }

        String getLabel() {
            return tokens[0];
        }

        String getType() {
            return tokens[1];
        }

        String getDefaultValue() {
            return tokens[2];
        }

        String getComment() {
            return tokens[3];
        }

        String getGroup() {
            return group;
        }

    }


    static LinkedHashMap<String, List<Preference>> loadPreferenceList() throws IOException {

        LinkedHashMap<String, List<Preference>> map = new LinkedHashMap<>();
        List<Preference> prefList = null;

        BufferedReader reader = new BufferedReader(new java.io.InputStreamReader(PreferenceEditorFX.class.getResourceAsStream("/resources/preferences.txt")));
        String nextLine;
        String group = null;
        Preference currentPreference = null;

        while ((nextLine = reader.readLine()) != null) {

            if (nextLine.startsWith("//")) continue;

            if (nextLine.startsWith("#Hidden")) break;

            nextLine = nextLine.trim();

            if (nextLine.startsWith(SEPARATOR_KEY)) {
                prefList.add(new Preference(SEPARATOR_KEY, group));   // "Blank" preference
                continue;
            }

            if (nextLine.startsWith("##")) {

                // New group

                if (prefList != null && currentPreference != null) {
                    prefList.add(currentPreference);
                }

                group = null;  // End previous group
                if (nextLine.length() > 2) {
                    group = nextLine.substring(2);
                }

                currentPreference = null;
                continue;
            }

            if (nextLine.startsWith("#")) {

                // New tab

                group = null;

                if (prefList != null && currentPreference != null) {
                    prefList.add(currentPreference);
                }

                String tab = nextLine.substring(1);
                prefList = new ArrayList<>();
                map.put(tab, prefList);
                currentPreference = null;
                continue;
            }

            if (nextLine.length() == 0) {
                if (prefList != null && currentPreference != null) {
                    prefList.add(currentPreference);
                    currentPreference = null;
                }
            } else if (currentPreference == null) {
                currentPreference = new Preference(nextLine, group);
            } else {
                currentPreference.pushToken(nextLine);
            }


        }
        return map;
    }

}
