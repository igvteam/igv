package org.broad.igv.prefs;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.StackPane;
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

        for(Map.Entry<String, List<Preference>> entry : prefMap.entrySet()) {

            Tab tab = new Tab(entry.getKey());
            tab.setClosable(false);

            pane.getTabs().add(tab);

            ScrollPane scrollPane = new ScrollPane();
            tab.setContent(scrollPane);

            GridPane gridPane = new GridPane();

            gridPane.setHgap(5);
            gridPane.setVgap(5);
            scrollPane.setContent(gridPane);

            int row = 1;
            for (Preference pref : entry.getValue()) {

                if (pref.type.equals("boolean")) {
                    CheckBox cb = new CheckBox(pref.label);
                    cb.setSelected(preferences.getAsBoolean(pref.key));
                    cb.setOnAction(event -> {
                        updatedPrefs.put(pref.key, Boolean.toString(cb.isSelected()));
                        System.out.println("Set " + pref.label + ": " + cb.isSelected());
                    });
                    GridPane.setColumnSpan(cb, 2);
                    gridPane.add(cb, 1, row);
                } else {
                    Label label = new Label(pref.label);
                    TextField field = new TextField(preferences.get(pref.key));
                    field.setOnAction(event -> {
                        final String text = field.getText();
                        if (validate(text, pref.type)) {
                            System.out.println("Set " + pref.label + ": " + text);
                            updatedPrefs.put(pref.key, text);
                        } else {
                            field.setText(preferences.get(pref.key));
                        }

                    });
                    gridPane.add(label, 1, row);
                    gridPane.add(field, 2, row);
                }

                if (pref.comment != null) {
                    Text label = new Text(pref.comment);
                    label.setStyle("-fx-font-style: italic;");  // Currently doesn't work (JavaFX bug)
                    gridPane.add(label, 3, row);
                }

                row++;
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

        String key;
        String label;
        String type;
        String comment;

        public Preference(String key, String label, String type, String comment) {
            this.comment = comment;
            this.key = key;
            this.label = label;
            this.type = type;
        }
    }


    static LinkedHashMap<String, List<Preference>> loadPreferenceList() throws IOException {

        LinkedHashMap<String, List<Preference>> map = new LinkedHashMap<>();
        List<Preference> prefList = null;

        BufferedReader reader = new BufferedReader(new java.io.InputStreamReader(PreferenceEditorFX.class.getResourceAsStream("/resources/preferences.txt")));
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            nextLine = nextLine.trim();

            if (nextLine.startsWith("//")) continue;

            if (nextLine.startsWith("#")) {
                String groupName = nextLine.substring(1);
                prefList = new ArrayList<>();
                map.put(groupName, prefList);

            } else if(prefList != null)  {
                String[] tokens = Globals.tabPattern.split(nextLine);

                if (tokens.length >= 3) {
                    String key = tokens[0];
                    String label = tokens[1];
                    String type = tokens[2].toLowerCase();
                    String comment = tokens.length > 3 ? tokens[3] : null;
                    prefList.add(new Preference(key, label, type, comment));
                }
            }
        }
        return map;
    }

}
