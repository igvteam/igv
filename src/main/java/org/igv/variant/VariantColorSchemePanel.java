package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.IGV;
import org.igv.ui.util.FileDialogUtils;
import org.igv.ui.util.MessageUtils;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.table.AbstractTableModel;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.io.File;
import java.util.List;

/**
 * Manages the color schemes used to color variants by a VCF INFO attribute.  Presented on the "Variants" tab of
 * the preferences editor.
 * <p>
 * The set of INFO keys, and of values for each key, is unbounded, so there is nothing to enumerate as individual
 * color preferences.  Instead schemes are files, and this panel manages which ones IGV has.  Importing copies the
 * file into the IGV directory, which is then the only copy IGV reads -- editing or moving the original has no
 * further effect.
 */
public class VariantColorSchemePanel extends JPanel {

    private static Logger log = LogManager.getLogger(VariantColorSchemePanel.class);

    private final SchemeTableModel tableModel = new SchemeTableModel();
    private final JTable table = new JTable(tableModel);
    private final JButton removeButton = new JButton("Remove");

    public VariantColorSchemePanel() {

        setLayout(new BorderLayout(0, 5));
        setBorder(BorderFactory.createTitledBorder("Variant INFO colors"));

        JLabel description = new JLabel("<html>Colors for the values of VCF INFO attributes, used by the track menu's "
                + "<i>Color By &gt; INFO Field</i> option.<br>Imported schemes are copied to "
                + VariantColorSchemes.getSchemeDirectory().getAbsolutePath() + ".");
        add(description, BorderLayout.NORTH);

        table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        table.setFillsViewportHeight(true);
        table.getSelectionModel().addListSelectionListener(e -> updateButtonState());
        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setPreferredSize(new Dimension(650, 110));
        add(scrollPane, BorderLayout.CENTER);

        JButton importFileButton = new JButton("Import File...");
        importFileButton.addActionListener(e -> importFile());

        JButton importUrlButton = new JButton("Import URL...");
        importUrlButton.addActionListener(e -> importUrl());

        removeButton.addActionListener(e -> removeSelected());

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        buttonPanel.add(importFileButton);
        buttonPanel.add(importUrlButton);
        buttonPanel.add(removeButton);
        add(buttonPanel, BorderLayout.SOUTH);

        updateButtonState();
    }

    private void importFile() {
        File file = FileDialogUtils.chooseFile("Select variant color scheme");
        if (file == null) {
            return;
        }
        try {
            VariantColorSchemes.importFile(file);
            schemesChanged();
        } catch (Exception e) {
            log.error("Error importing variant color scheme: " + file.getAbsolutePath(), e);
            MessageUtils.showMessage("Error importing color scheme: " + e.getMessage());
        }
    }

    private void importUrl() {
        String url = MessageUtils.showInputDialog("Color scheme URL");
        if (url == null || url.trim().isEmpty()) {
            return;
        }
        try {
            VariantColorSchemes.importUrl(url.trim());
            schemesChanged();
        } catch (Exception e) {
            log.error("Error importing variant color scheme: " + url, e);
            MessageUtils.showMessage("Error importing color scheme: " + e.getMessage());
        }
    }

    private void removeSelected() {
        VariantColorScheme scheme = getSelectedScheme();
        if (scheme == null || scheme.isBuiltIn()) {
            return;
        }
        if (MessageUtils.confirm("Remove color scheme \"" + scheme.getName() + "\"?  Its file will be deleted.")) {
            VariantColorSchemes.remove(scheme);
            schemesChanged();
        }
    }

    private void schemesChanged() {
        tableModel.refresh();
        updateButtonState();
        // Colors resolve through the schemes on every repaint, so open tracks pick this up immediately
        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
    }

    private void updateButtonState() {
        VariantColorScheme selected = getSelectedScheme();
        removeButton.setEnabled(selected != null && !selected.isBuiltIn());
    }

    private VariantColorScheme getSelectedScheme() {
        int row = table.getSelectedRow();
        return row < 0 ? null : tableModel.getScheme(table.convertRowIndexToModel(row));
    }

    private static class SchemeTableModel extends AbstractTableModel {

        private static final String[] COLUMNS = {"Name", "INFO attributes", "Source"};

        private List<VariantColorScheme> schemes = VariantColorSchemes.getSchemes();

        void refresh() {
            schemes = VariantColorSchemes.getSchemes();
            fireTableDataChanged();
        }

        VariantColorScheme getScheme(int row) {
            return row >= 0 && row < schemes.size() ? schemes.get(row) : null;
        }

        @Override
        public int getRowCount() {
            return schemes.size();
        }

        @Override
        public int getColumnCount() {
            return COLUMNS.length;
        }

        @Override
        public String getColumnName(int column) {
            return COLUMNS[column];
        }

        @Override
        public Object getValueAt(int row, int column) {
            VariantColorScheme scheme = schemes.get(row);
            return switch (column) {
                case 0 -> scheme.getName();
                case 1 -> String.join(", ", scheme.getKeys());
                case 2 -> scheme.isBuiltIn() ? "Built in" :
                        (scheme.getSource() == null ? scheme.getFile().getName() : scheme.getSource());
                default -> "";
            };
        }
    }
}
