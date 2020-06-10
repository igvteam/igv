package org.broad.igv.prefs;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.google.OAuthUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorSwatch;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.Utilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

public class PreferencesEditor {

    private static Logger log = Logger.getLogger(PreferencesEditor.class);

    private static final Font labelFont = new Font("Lucida Grande", Font.BOLD, 14);

    public static void main(String[] args) throws IOException {
        open(null);
    }

    public static void open(Frame parent) throws IOException {
        List<PreferencesManager.PreferenceGroup> preferenceGroups = PreferencesManager.loadPreferenceList();

        SwingUtilities.invokeLater(() -> {
            JDialog frame = new JDialog(parent, "Preferences", true);
            final JPanel panel = new JPanel();

            SwingUtilities.invokeLater(() -> {
                init(frame, panel, preferenceGroups);
            });
            panel.setPreferredSize(new Dimension(850, 590));
            panel.setMaximumSize(new Dimension(850, 590));
            frame.add(panel);
            frame.pack();
            frame.setSize(900, 600);
            frame.setLocationRelativeTo(parent);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        });
    }

    private static void init(final JDialog parent, JPanel panel, List<PreferencesManager.PreferenceGroup> preferenceGroups) {

        // final Map<String, String> updatedPrefs = new HashMap<>();

        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.setPreferredSize(new Dimension(850, 590));
        tabbedPane.setMaximumSize(new Dimension(850, 590));
        panel.setLayout(new BorderLayout());
        panel.add(tabbedPane, BorderLayout.CENTER);

        Map<String, Map<String, String>> updatedPreferencesMap = new HashMap<>();
        for (PreferencesManager.PreferenceGroup entry : preferenceGroups) {

            if (entry.tabLabel.equals("Hidden")) continue;

            final IGVPreferences preferences = PreferencesManager.getPreferences(entry.category);

            final Map<String, String> updatedPrefs =
                    updatedPreferencesMap.containsKey(entry.category) ? updatedPreferencesMap.get(entry.category) :
                            new HashMap<>();
            updatedPreferencesMap.put(entry.category, updatedPrefs);

            final String tabLabel = entry.tabLabel;
            JScrollPane scrollPane = new JScrollPane();
            scrollPane.setPreferredSize(new Dimension(750, 590));
            scrollPane.setMaximumSize(new Dimension(750, 590));
            scrollPane.setName(tabLabel);
            scrollPane.getVerticalScrollBar().setUnitIncrement(16);
            tabbedPane.addTab(tabLabel, scrollPane);

            JPanel content = new JPanel();
//            content.setLayout(new BoxLayout(content, BoxLayout.Y_AXIS));
            int contentRow = 0;
            GridBagLayout contentGrid = new GridBagLayout();
            content.setLayout(contentGrid);
            scrollPane.setViewportView(content);

            JPanel group = new JPanel();

            int row = 0;
            GridBagLayout grid = new GridBagLayout();
            group.setLayout(grid);
            contentGrid.addLayoutComponent(group, new GridBagConstraints(0, contentRow++, 1, 1, 1.0, 0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 5, 5));
            content.add(group);

            String currentGroup = null;

            // Add an empty spacer component to fill any remaining horizontal group space
            JLabel groupSpacer = new JLabel("    ");
            groupSpacer.setPreferredSize(new Dimension(1, 1));
            grid.addLayoutComponent(groupSpacer, new GridBagConstraints(0, row++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
            group.add(groupSpacer);

            for (PreferencesManager.Preference pref : entry.preferences) {

                try {
                    if (pref.getKey().equals(PreferencesManager.SEPARATOR_KEY)) {
                        JSeparator sep = new JSeparator();
                        grid.addLayoutComponent(sep, new GridBagConstraints(0, row, 4, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 5, 2));
                        group.add(sep);
                        row++;
                        continue;
                    }

                    if (pref.getKey().equals(PreferencesManager.INFO_KEY)) {
                        row++;
                        JLabel label = new JLabel(pref.getLabel());
                        label.setFont(labelFont);
                        grid.addLayoutComponent(label, new GridBagConstraints(0, row, 4, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(3, 5, 5, 5), 5, 2));
                        group.add(label);
                        row++;
                        continue;
                    }

                    if (pref.group != null && !pref.group.equals(currentGroup)) {
                        // Start a new group
                        row = 0;

                        if (group.getComponentCount() > 0) {
                            // Don't create a new JPanel and Layout if nothing has been added to the existing one.
                            group = new JPanel();
                            grid = new GridBagLayout();
                            group.setLayout(grid);
                            contentGrid.addLayoutComponent(group, new GridBagConstraints(0, contentRow++, 1, 1, 1.0, 0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(10, 10, 10, 10), 5, 5));
                            content.add(group);

                            // Add an empty spacer component to fill any remaining horizontal group space
                            groupSpacer = new JLabel("    ");
                            groupSpacer.setPreferredSize(new Dimension(1, 1));
                            grid.addLayoutComponent(groupSpacer, new GridBagConstraints(0, row++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
                            group.add(groupSpacer);
                        }

                        currentGroup = pref.group;
                        group.setBorder(BorderFactory.createTitledBorder(currentGroup));
                    }

                    if (pref.getType().equals("boolean")) {
                        JCheckBox cb = new JCheckBox(pref.getLabel());
                        cb.setSelected(preferences.getAsBoolean(pref.getKey()));
                        cb.addActionListener(event -> {
                            updatedPrefs.put(pref.getKey(), Boolean.toString(cb.isSelected()));
                        });

                        grid.addLayoutComponent(cb, new GridBagConstraints(0, row, 2, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(3, 5, 2, 5), 2, 2));
                        group.add(cb);

                        if (pref.getComment() != null) {
                            cb.setToolTipText(pref.getComment());
                        }

                    } else if (pref.getType().startsWith("select")) {
                        JLabel label = new JLabel(pref.getLabel());
                        String[] selections = Globals.whitespacePattern.split(pref.getType())[1].split("\\|");
                        final JComboBox<String> comboBox = new JComboBox<String>(selections);
                        comboBox.setSelectedItem(pref.getDefaultValue().toString());
                        comboBox.addActionListener(event -> {
                            log.debug("Set " + pref.getLabel() + " " + comboBox.getSelectedItem());
                        });
                        grid.addLayoutComponent(label, new GridBagConstraints(0, row, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 5, 2, 3), 2, 2));
                        grid.addLayoutComponent(comboBox, new GridBagConstraints(1, row, 3, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 2, 2, 5), 2, 2));
                        group.add(label);
                        group.add(comboBox);

                        if (pref.getComment() != null) {
                            label.setToolTipText(pref.getComment());
                            comboBox.setToolTipText(pref.getComment());
                        }
                    } else if(pref.getType().startsWith("color")) {
                        String colorString = preferences.get(pref.getKey());
                        Color c;
                        if(colorString == null) {
                            c = Color.WHITE;
                        } else {
                            c = ColorUtilities.stringToColor(colorString);
                        }
                        JLabel label = new JLabel(pref.getLabel());
                        ColorSwatch colorSwatch = new ColorSwatch(c);
                        colorSwatch.addColorChangeListener(c1 -> {
                            String s = ColorUtilities.colorToString(c1);
                            updatedPrefs.put(pref.getKey(), s);
                        });
                        grid.addLayoutComponent(label, new GridBagConstraints(0, row, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 5, 2, 3), 2, 2));
                        grid.addLayoutComponent(colorSwatch, new GridBagConstraints(1, row, 3, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 2, 2, 5), 2, 2));
                        group.add(label);
                        group.add(colorSwatch);

                    }

                    else {
                        JLabel label = new JLabel(pref.getLabel());
                        JTextField field = new JTextField(preferences.get(pref.getKey()));
                        Dimension d = field.getPreferredSize();
                        d.width = 300;
                        field.setPreferredSize(d);
                        field.setMaximumSize(d);
                        field.addActionListener(event -> {
                            String text = field.getText();
                            if (validate(text, pref.getType())) {
                                // TODO -- make base64 an explicit type
                                if (pref.getKey().equals(Constants.PROXY_PW)) {
                                    text = Utilities.base64Encode(text);
                                }

                                updatedPrefs.put(pref.getKey(), text);
                            } else {
                                field.setText(preferences.get(pref.getKey()));
                            }
                        });
                        field.addFocusListener(new FocusAdapter() {
                            @Override
                            public void focusLost(FocusEvent e) {
                                // Validate and save the value if the Preference field loses focus
                                String text = field.getText();
                                if (validate(text, pref.getType())) {
                                    // TODO -- make base64 an explicit type
                                    if (pref.getKey().equals(Constants.PROXY_PW)) {
                                        text = Utilities.base64Encode(text);
                                    }
                                    updatedPrefs.put(pref.getKey(), text);
                                } else {
                                    field.setText(preferences.get(pref.getKey()));
                                }
                            }
                        });

                        grid.addLayoutComponent(label, new GridBagConstraints(0, row, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 5, 2, 3), 2, 2));
                        grid.addLayoutComponent(field, new GridBagConstraints(1, row, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 2, 5), 2, 2));
                        group.add(label);
                        group.add(field);

                        if (pref.getComment() != null) {
                            label.setToolTipText(pref.getComment());
                            field.setToolTipText(pref.getComment());
                        }
                    }

                    row++;
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            if (tabLabel.equalsIgnoreCase("Cram")) {
                // Add Cram cache directory management at the end.  This is a special case

                String currentDirectory = DirectoryManager.getFastaCacheDirectory().getAbsolutePath();
                final JLabel currentDirectoryLabel = new JLabel("Cache directory: " + currentDirectory);
                final JButton moveButton = new JButton("Move...");
                row++;
                grid.addLayoutComponent(currentDirectoryLabel, new GridBagConstraints(0, row, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 5, 2, 5), 2, 2));
                grid.addLayoutComponent(moveButton, new GridBagConstraints(1, row, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 5, 2, 5), 2, 2));
                group.add(currentDirectoryLabel);
                group.add(moveButton);

                moveButton.addActionListener(event -> {
                    UIUtilities.invokeOnEventThread(() -> {
                        final File directory = DirectoryManager.getFastaCacheDirectory();
                        final File newDirectory = FileDialogUtils.chooseDirectory("Select cache directory", DirectoryManager.getUserDirectory());
                        if (newDirectory != null && !newDirectory.equals(directory)) {
                            DirectoryManager.moveDirectoryContents(directory, newDirectory);
                            SwingUtilities.invokeLater(() -> currentDirectoryLabel.setText(newDirectory.getAbsolutePath()));
                        }
                    });
                });

            }


            if (tabLabel.equalsIgnoreCase("Advanced")) {
                // Add IGV directory management at the end.  This is a special case
                String currentDirectory = DirectoryManager.getIgvDirectory().getAbsolutePath();
                final JLabel currentDirectoryLabel = new JLabel("IGV Directory: " + currentDirectory);
                final JButton moveButton = new JButton("Move...");
                row++;
                grid.addLayoutComponent(currentDirectoryLabel, new GridBagConstraints(0, row, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(3, 5, 2, 5), 2, 2));
                grid.addLayoutComponent(moveButton, new GridBagConstraints(1, row, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 5, 2, 5), 2, 2));
                group.add(currentDirectoryLabel);
                group.add(moveButton);

                moveButton.addActionListener(event -> {
                    UIUtilities.invokeOnEventThread(() -> {
                        final File igvDirectory = DirectoryManager.getIgvDirectory();
                        final File newDirectory = FileDialogUtils.chooseDirectory("Select IGV directory", DirectoryManager.getUserDirectory());
                        if (newDirectory != null && !newDirectory.equals(igvDirectory)) {
                            DirectoryManager.moveIGVDirectory(newDirectory);
                            SwingUtilities.invokeLater(() -> currentDirectoryLabel.setText(newDirectory.getAbsolutePath()));
                        }
                    });
                });

            }

            // Add an empty spacer component to fill any remaining vertical content space in order to give a compact presentation.  This will occupy no appreciable space unless contents are too small to fill the ScrollPane
            JPanel contentSpacer = new JPanel();
            contentGrid.addLayoutComponent(contentSpacer, new GridBagConstraints(0, contentRow++, 1, 1, 1.0, 1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 2, 2));
            content.add(contentSpacer);
        }


        JPanel saveCancelPanel = new JPanel();
        saveCancelPanel.setBackground(new Color(0x33, 0x66, 0x99));
        saveCancelPanel.setPreferredSize(new Dimension(750, 40));
        saveCancelPanel.setMaximumSize(new Dimension(750, 40));

        JButton cancelButton = new JButton("Cancel");
        cancelButton.setPreferredSize(new Dimension(100, 30));
        cancelButton.setMaximumSize(new Dimension(100, 30));
        cancelButton.addActionListener((event) -> {
            SwingUtilities.invokeLater(() -> parent.setVisible(false));
        });

        JButton saveButton = new JButton("Save");
        saveButton.setPreferredSize(new Dimension(100, 30));
        saveButton.setMaximumSize(new Dimension(100, 30));
        saveButton.setDefaultCapable(true);
        saveButton.addActionListener((event) -> {
            saveAction(event, updatedPreferencesMap);
            SwingUtilities.invokeLater(() -> parent.setVisible(false));
        });
        saveCancelPanel.add(cancelButton);
        saveCancelPanel.add(saveButton);

        panel.add(saveCancelPanel, BorderLayout.SOUTH);
    }

    private static void saveAction(ActionEvent event, Map<String, Map<String, String>> updatedPreferencesMap) {

        for (Map<String, String> prefs : updatedPreferencesMap.values()) {
            extractMutationPreferences(prefs);
        }

        PreferencesManager.updateAll(updatedPreferencesMap);

        for (Map<String, String> map : updatedPreferencesMap.values()) {
            if (map.containsKey(PROVISIONING_URL)) {
                try {
                    OAuthUtils.getInstance().loadProvisioningURL(map.get(PROVISIONING_URL));
                } catch (IOException e) {
                    MessageUtils.showErrorMessage("Error loading provisioning URL", e);
                }
                break;
            }
        }

        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
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


    static private void extractMutationPreferences(Map<String, String> prefs) {

        PaletteColorTable ct = PreferencesManager.getPreferences().getMutationColorScheme();
        Set<String> keys = new HashSet(ct.getKeys());
        keys.addAll(Arrays.asList("Indel", "Missense", "Nonsense", "Splice_site", "Synonymous", "Targeted_Region",
                "Unknown", "Truncating", "Non-coding_Transcript", "Other_AA_changing", "Other_likely_neutral"));

        List<Map.Entry<String, String>> mutColors = new ArrayList<>();
        for (Map.Entry<String, String> entry : prefs.entrySet()) {
            if (keys.contains(entry.getKey())) {
                mutColors.add(entry);
            }
        }

        if (!mutColors.isEmpty()) {
            for (Map.Entry<String, String> entry : mutColors) {
                ct.getColorMap().put(entry.getKey(), ColorUtilities.stringToColor(entry.getValue()));
            }
            String mapString = ct.getMapAsString();
            prefs.put(MUTATION_COLOR_TABLE, mapString);

            for (String key : keys) {
                prefs.remove(key);
            }
        }
    }
}
