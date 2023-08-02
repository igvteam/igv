package org.broad.igv.ui.form;

import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesTextField;
import org.broad.igv.ui.color.ColorSwatch;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.Utilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.broad.igv.ui.color.ColorUtilities.*;

public class FormEditor {

    private static Logger log = LogManager.getLogger(FormEditor.class);

    private static final Font labelFont = new Font("Lucida Grande", Font.BOLD, 14);
    public static final String SEPARATOR_KEY = "---";
    public static final String INFO_KEY = "info";


    public static void main(String[] args) throws IOException {

        Map<String, java.util.List<FormField>> preferences = new HashMap<>();

        preferences.put("", Arrays.asList(
                new FormField("threshold", "Threshold", "float", "0.5", ""),
                new FormField("twocolor", "Two Color", "boolean", "false", "Color scheme")
        ));
        preferences.put("Show modifications", Arrays.asList(
                new FormField("all", "All", "boolean", "true", "")

        ));


        open(null, preferences);
    }

    public static void open(Frame parent, Map<String, java.util.List<FormField>> preferences) throws IOException {


        SwingUtilities.invokeLater(() -> {
            JDialog dialog = new JDialog(parent, "Preferences", true);

            SwingUtilities.invokeLater(() -> {
                init(dialog, preferences);
            });
            //panel.setPreferredSize(new Dimension(850, 590));
            //panel.setMaximumSize(new Dimension(850, 590));

            dialog.setPreferredSize(new Dimension(850, 590));
            dialog.pack();
            //frame.setSize(900, 600);
            dialog.setLocationRelativeTo(parent);
            dialog.setVisible(true);
            dialog.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        });
    }

    private static void init(final JDialog dialog, Map<String, java.util.List<FormField>> preferences) {

        Map<String, String> updatedPrefs = new HashMap<>();

        //panel.setLayout(new BorderLayout());

        JPanel content = new JPanel();
        content.setLayout(new BorderLayout());
        dialog.add(content);

        JPanel formContent = new JPanel();
        formContent.setLayout(new BoxLayout(formContent, BoxLayout.PAGE_AXIS));
        content.add(formContent, BorderLayout.CENTER);

        int row = 0;
        for (Map.Entry<String, java.util.List<FormField>> entry : preferences.entrySet()) {

            String currentGroup = entry.getKey();
            // Start a new group
            row = 0;

            JPanel group = new JPanel();
            group.setLayout(new BoxLayout(group, BoxLayout.PAGE_AXIS));
            formContent.add(group);

            group.setBorder(BorderFactory.createTitledBorder(currentGroup));

            for (FormField pref : entry.getValue()) {

                JPanel fieldPanel = new JPanel();
                fieldPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
                fieldPanel.setMaximumSize(new Dimension(1000, 30));
                group.add(fieldPanel);

                if (pref.getKey().equals(SEPARATOR_KEY)) {
                    JSeparator sep = new JSeparator();
                    group.add(sep);
                    row++;
                    continue;
                }

                if (pref.getKey().equals(INFO_KEY)) {
                    row++;
                    JLabel label = new JLabel(pref.getLabel());
                    label.setFont(labelFont);
                    group.add(label);
                    row++;
                    continue;
                }


                if (pref.getType().equals("boolean")) {
                    JCheckBox cb = new JCheckBox(pref.getLabel());
//                        cb.setSelected(preferences.getAsBoolean(pref.getKey()));
//                        cb.addActionListener(event -> {
//                            updatedPrefs.put(pref.getKey(), Boolean.toString(cb.isSelected()));
//                        });


                    if (pref.getComment() != null) {
                        cb.setToolTipText(pref.getComment());
                    }

                    fieldPanel.add(cb);

                } else if (pref.getType().startsWith("select")) {
                    JLabel label = new JLabel(pref.getLabel());

                    String[] selections = Globals.whitespacePattern.split(pref.getType())[1].split("\\|");
//                        String[] labels = Arrays.stream(selections)
//                                .map((s) -> labelMappings.containsKey(s) ? labelMappings.get(s) : s)
//                                .toArray(String[]::new);


                    final JComboBox<String> comboBox = new JComboBox<String>(selections);
//                        comboBox.setSelectedItem(preferences.get(pref.getKey()));
//                        comboBox.addActionListener(event -> {
//                            updatedPrefs.put(pref.getKey(), comboBox.getSelectedItem().toString());
//                        });

                    fieldPanel.add(label);
                    fieldPanel.add(comboBox);

                    if (pref.getComment() != null) {
                        label.setToolTipText(pref.getComment());
                        comboBox.setToolTipText(pref.getComment());
                    }
                } else if (pref.getType().startsWith("color")) {
                    String colorString = pref.getDefaultValue();
                    Color c;
                    if (colorString == null) {
                        c = Color.WHITE;
                    } else {
                        c = stringToColor(colorString);
                    }
                    JLabel label = new JLabel(pref.getLabel());
                    ColorSwatch colorSwatch = new ColorSwatch(c);
                    colorSwatch.addColorChangeListener(c1 -> {
                        String s = colorToString(c1);
                        updatedPrefs.put(pref.getKey(), s);
                    });

                    fieldPanel.add(label);
                    fieldPanel.add(colorSwatch);

                } else {
                    JLabel label = new JLabel(pref.getLabel());

                    String fieldText = pref.getDefaultValue();
                    if (pref.getKey().equals(Constants.PROXY_PW) && fieldText != null && fieldText.length() > 0) {
                        fieldText = Utilities.base64Decode(fieldText);
                    }
                    PreferencesTextField field = new PreferencesTextField(pref.getKey().equals(Constants.PROXY_PW) ? new JPasswordField(fieldText) : new JTextField(fieldText));
                    Dimension d = field.get().getPreferredSize();
                    d.width = 300;
                    field.get().setPreferredSize(d);
                    field.get().setMaximumSize(d);
                    field.get().addActionListener(event -> {
                        String text = field.getPreferenceText();
                        if (validate(text, pref.getType())) {
                            if (pref.getKey().equals(Constants.PROXY_PW)) {
                                text = Utilities.base64Encode(text);
                            }
                            updatedPrefs.put(pref.getKey(), text);
                        } else {
                            MessageUtils.showMessage("Field " + pref.getLabel() + " must be a " + pref.getType());
                            field.get().setText(pref.getDefaultValue());
                        }
                    });
                    field.get().addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            // Validate and save the value if the Preference field loses focus
                            String text = field.getPreferenceText();
                            if (validate(text, pref.getType())) {
                                // TODO -- make base64 an explicit type
                                if (pref.getKey().equals(Constants.PROXY_PW)) {
                                    text = Utilities.base64Encode(text);
                                }
                                updatedPrefs.put(pref.getKey(), text);
                            } else {
                                MessageUtils.showMessage("Field " + pref.getLabel() + " must be a " + pref.getType());
                                field.get().setText(pref.getDefaultValue());
                            }
                        }
                    });


                    fieldPanel.add(label);
                    fieldPanel.add(field.get());

                    if (pref.getComment() != null) {
                        label.setToolTipText(pref.getComment());
                        field.get().setToolTipText(pref.getComment());
                    }
                }

                row++;

            }
        }


        JPanel saveCancelPanel = new JPanel();
        saveCancelPanel.setBackground(new Color(0x33, 0x66, 0x99));
        saveCancelPanel.setPreferredSize(new Dimension(750, 40));
        saveCancelPanel.setMaximumSize(new Dimension(750, 40));

        JButton cancelButton = new JButton("Cancel");
        cancelButton.setPreferredSize(new Dimension(100, 30));
        cancelButton.setMaximumSize(new Dimension(100, 30));
        cancelButton.addActionListener((event) -> {
            SwingUtilities.invokeLater(() -> dialog.setVisible(false));
        });

        JButton saveButton = new JButton("Save");
        saveButton.setPreferredSize(new Dimension(100, 30));
        saveButton.setMaximumSize(new Dimension(100, 30));
        saveButton.setDefaultCapable(true);
        saveButton.addActionListener((event) -> {
            saveAction(event, updatedPrefs);
            SwingUtilities.invokeLater(() -> dialog.setVisible(false));
        });
        saveCancelPanel.add(cancelButton);
        saveCancelPanel.add(saveButton);

        content.add(saveCancelPanel, BorderLayout.SOUTH);
    }

    private static void saveAction(ActionEvent event, Map<String, String> updatedPreferencesMap) {

//        for (Map<String, String> prefs : updatedPreferencesMap.values()) {
//            extractMutationPreferences(prefs);
//        }
//
//        PreferencesManager.updateAll(updatedPreferencesMap);
//
//        for (Map<String, String> map : updatedPreferencesMap.values()) {
//            if (map.containsKey(PROVISIONING_URL)) {
//                try {
//                    OAuthUtils.getInstance().updateOauthProvider(map.get(PROVISIONING_URL));
//                } catch (IOException e) {
//                    MessageUtils.showErrorMessage("Error loading provisioning URL", e);
//                }
//                break;
//            }
//        }
//
//        if (IGV.hasInstance()) {
//            IGV.getInstance().repaint();
//        }
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
