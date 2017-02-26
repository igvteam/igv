package org.broad.igv.prefs;

import javafx.application.Application;
import javafx.application.Platform;
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
import org.broad.igv.Globals;

import javax.swing.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

public class PreferenceEditorFX {

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

            Tab tab = new Tab(entry.getKey());
            tab.setClosable(false);

            pane.getTabs().add(tab);

            ScrollPane scrollPane = new ScrollPane();
            tab.setContent(scrollPane);

            VBox vBox = new VBox();
            vBox.setFillWidth(true);
            vBox.prefWidthProperty().bind(scene.widthProperty());
            scrollPane.setContent(vBox);

            GridPane gridPane = new GridPane();
            gridPane.setHgap(5);
            gridPane.setVgap(5);
            vBox.getChildren().add(gridPane);

            String currentGroup = null;

            int row = 1;
            for (Preference pref : entry.getValue()) {

                try {
                    if (pref.getKey() == null) {
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
        String[] tokens;


        Preference(String group) {
            tokens = new String[6];
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
            return tokens[0];
        }

        String getLabel() {
            return tokens[1];
        }

        String getType() {
            return tokens[2];
        }

        String getDefaultValue() {
            return tokens[3];
        }

        String getComment() {
            return tokens[4];
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

            if (nextLine.startsWith("---")) {
                prefList.add(new Preference(group));   // "Blank" preference
                continue;
            }

            if (nextLine.startsWith("##")) {
                if(prefList != null && currentPreference != null) {
                    prefList.add(currentPreference);
                }

                group = null;  // End previous group
                if (nextLine.length() > 2) {
                    group = nextLine.substring(2);
                }

                currentPreference = new Preference(group);
                continue;
            }

            if (nextLine.startsWith("#")) {

                group = null;

                if(prefList != null && currentPreference != null) {
                    prefList.add(currentPreference);
                }

                String tab = nextLine.substring(1);
                prefList = new ArrayList<>();
                map.put(tab, prefList);
                currentPreference = new Preference(group);
                continue;
            }

            if(nextLine.length() == 0) {
                if(prefList != null && currentPreference != null) {
                    prefList.add(currentPreference);
                }
                currentPreference = new Preference(group);
            } else if(currentPreference != null) {
                currentPreference.pushToken(nextLine);
            }


        }
        return map;
    }

}
