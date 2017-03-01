package org.broad.igv.prefs;

import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.embed.swing.JFXPanel;
import javafx.geometry.HPos;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.text.Text;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class PreferenceEditorFX {

    public static final String SEPARATOR_KEY = "---";
    public static final String INFO_KEY = "info";
    public static final String CATEGORY_KEY = "category";

    public static void main(String[] args) throws IOException {

        List<PreferenceGroup> preferenceGroups = loadPreferenceList();

        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Preferences");
            final JFXPanel fxPanel = new JFXPanel();
            frame.add(fxPanel);
            frame.setSize(800, 600);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

            Platform.runLater(() -> initFX(fxPanel, preferenceGroups));
        });
    }

    private static void initFX(JFXPanel fxPanel, List<PreferenceGroup> preferenceGroups) {

        final Map<String, String> updatedPrefs = new HashMap<>();

        TabPane pane = new TabPane();
        Scene scene = new Scene(pane);

        for (PreferenceGroup entry : preferenceGroups) {

            final IGVPreferences preferences = PreferencesManager.getPreferences(entry.category);

            final String tabLabel = entry.tabLabel;
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
            for (Preference pref : entry.preferences) {

                try {
                    Tooltip tooltip = pref.getComment() == null ? null : new Tooltip(pref.getComment());

                    if (pref.getKey().equals(SEPARATOR_KEY)) {
                        Separator sep = new Separator();
                        GridPane.setColumnSpan(sep, 4);
                        gridPane.add(sep, 1, row);
                        row++;
                        continue;
                    }

                    if (pref.getKey().equals(INFO_KEY)) {
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
                                System.out.println("Set " + pref.getLabel() + ": " + text);
                                updatedPrefs.put(pref.getKey(), text);
                            } else {
                                field.setText(preferences.get(pref.getKey()));
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

        String group;
        String[] tokens;

        Preference(String [] tokens, String group) {
            this.tokens = tokens;
            this.group = group;
        }

        Preference (String key, String group)
        {
            this(new String [] {key, null, null, null}, group);
        }

        Preference (String key, String label, String group)
        {
            this(new String [] {key, label, null, null}, group);
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
            return tokens.length > 4 ? tokens[4] : null;
        }

        String getGroup() {
            return group;
        }

        String printString() {
            String str = getKey() + "\t" + getLabel() + "\t" + getType() + "\t" + getDefaultValue();
            if (getComment() != null) str += "\t" + getComment();
            return str;
        }

    }


    static List<PreferenceGroup> loadPreferenceList() throws IOException {

        List<PreferenceGroup> groupList = new ArrayList<>();
        List<Preference> prefList = null;

        BufferedReader reader = new BufferedReader(new java.io.InputStreamReader(PreferenceEditorFX.class.getResourceAsStream("/org/broad/igv/prefs/preferences.tab")));
        String nextLine;
        String group = null;

        while ((nextLine = reader.readLine()) != null) {
            nextLine = nextLine.trim();

            if (nextLine.startsWith("//") || nextLine.length() == 0) {
                continue;
            }

            else if (nextLine.startsWith("#Hidden")) {
                break;
            }

            else if (nextLine.startsWith(SEPARATOR_KEY)) {
                prefList.add(new Preference(SEPARATOR_KEY, group));
                continue;
            }

            else if (nextLine.startsWith(INFO_KEY)) {
                Preference preference = new Preference(INFO_KEY, nextLine.substring(INFO_KEY.length()).trim(), group);
                prefList.add(preference);   // "Blank" preference
                continue;
            }


            else if (nextLine.startsWith("##")) {

                group = null;  // End previous group
                if (nextLine.length() > 2) {
                    group = nextLine.substring(2);  // New group
                }
                continue;
            }

            else if (nextLine.startsWith("#")) {

                // New tab

                String[] tokens = Globals.tabPattern.split(nextLine);
                String tabLabel = tokens[0].substring(1);
                String category = tokens.length > 1 ? tokens[1] : null;
                prefList = new ArrayList<>();
                PreferenceGroup preferenceGroup = new PreferenceGroup(tabLabel, category, prefList);
                groupList.add(preferenceGroup);

                group = null;

                continue;
            }

            else {

                String[] tokens = Globals.tabPattern.split(nextLine);
                if(tokens.length < 4) {

                }
                else {
                    prefList.add(new Preference(tokens, group));
                }
            }

        }

        return groupList;
    }

    static class PreferenceGroup {

        String tabLabel;
        String category;
        List<Preference> preferences;

        public PreferenceGroup(String tabLabel, String category, List<Preference> preferences) {
            this.tabLabel = tabLabel;
            this.category = category;
            this.preferences = preferences;
        }
    }

}
