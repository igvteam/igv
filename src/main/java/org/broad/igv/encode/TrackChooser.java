/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * Created by JFormDesigner on Thu Oct 31 22:31:02 EDT 2013
 */

package org.broad.igv.encode;

import com.jidesoft.swing.JideBoxLayout;
import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.Pair;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.TableColumn;
import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.text.ParseException;
import java.util.List;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * @author Jim Robinson
 */
public class TrackChooser extends org.broad.igv.ui.IGVDialog {

    private static Logger log = LogManager.getLogger(TrackChooser.class);

    private static NumberFormatter numberFormatter = new NumberFormatter();

    JTable table;
    JTextField filterTextField;
    JLabel rowCountLabel;
    TrackChooserModel model;
    private boolean canceled;


    public TrackChooser(Frame owner, final List<String> headings, final List<FileRecord> rows, String title) {
        super(owner);
        setTitle(title);
        this.model = new TrackChooserModel(headings, rows);
        setModal(true);
        initComponents(owner);
        init(model);
    }

    private void init(final TrackChooserModel model) {
        setModal(true);

        table.setAutoCreateRowSorter(true);
        table.setModel(model);
        table.setRowSorter(model.getSorter());
        try {
            rowCountLabel.setText(numberFormatter.valueToString(table.getRowCount()) + " rows");
        } catch (ParseException e) {
            log.error("Error parsing row count", e);
        }

        table.setRowSelectionAllowed(false);
        table.setColumnSelectionAllowed(false);

        TableColumn selectColumn = table.getColumnModel().getColumn(0);
        selectColumn.setPreferredWidth(25);
        selectColumn.setMaxWidth(30);

        filterTextField.getDocument().addDocumentListener(
                new DocumentListener() {
                    public void changedUpdate(DocumentEvent e) {
                        updateFilter();
                    }

                    public void insertUpdate(DocumentEvent e) {
                        updateFilter();
                    }

                    public void removeUpdate(DocumentEvent e) {
                        updateFilter();
                    }
                });

    }

    @Override
    public void setVisible(boolean b) {
        if (b) {
            this.model.updateSelections();
        }
        super.setVisible(b);
    }

    /**
     * Update the row filter regular expression from the expression in
     * the text box.
     */
    private void updateFilter() {


        RowFilter<TrackChooserModel, Object> rf = null;
        //If current expression doesn't parse, don't update.
        try {
            rf = new RegexFilter(filterTextField.getText());
        } catch (java.util.regex.PatternSyntaxException e) {
            return;
        }
        model.getSorter().setRowFilter(rf);

        try {
            rowCountLabel.setText(numberFormatter.valueToString(table.getRowCount()) + " rows");
        } catch (ParseException e) {
            log.error("Error parsing row count", e);
        }
    }


    private void loadButtonActionPerformed(ActionEvent e) {
        canceled = false;
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
    }

    public boolean isCanceled() {
        return canceled;
    }

    public List<FileRecord> getSelectedRecords()  {
        return model.getRecords().stream()
                .filter(record -> record.isSelected())
                .collect(Collectors.toUnmodifiableList());
    }

    public List<FileRecord> getAllRecords() {
        return model.getRecords();
    }


    private class RegexFilter extends RowFilter {

        List<Pair<String, Matcher>> matchers;

        RegexFilter(String text) {

            if (text == null) {
                throw new IllegalArgumentException("Pattern must be non-null");
            }
            matchers = new ArrayList<>();
            String[] tokens = Globals.whitespacePattern.split(text);
            for (String t : tokens) {
                // If token contains an = sign apply to specified column only
                String column = "*";
                String value = t.trim();
                if (t.contains("=")) {
                    String[] kv = Globals.equalPattern.split(t);
                    if (kv.length > 1) {
                        column = kv[0].trim();
                        value = kv[1].trim();
                    } else {
                        value = kv[0];  // Value is column name until more input is entered
                    }
                }

                matchers.add(new Pair(column, Pattern.compile("(?i)" + value).matcher("")));
            }

        }

        /**
         * Include row if each matcher succeeds in at least one column.  In other words all the conditions
         * are combined with "and"
         *
         * @param value
         * @return
         */
        @Override
        public boolean include(Entry value) {

            for (Pair<String, Matcher> entry : matchers) {
                String column = entry.getFirst();
                Matcher matcher = entry.getSecond();


                // Search for a match in at least one column.  The first column is the checkbox.
                boolean found = false;  // Pessimistic
                int nColumns = table.getColumnCount();
                for (int index = 1; index < nColumns; index++) {

                    // Include column headings in search.  This is to prevent premature filtering when entering a
                    // specific column condition (e.g. cataType=ChipSeq)
                    matcher.reset(table.getColumnName(index).toLowerCase());
                    if (matcher.find()) {
                        found = true;
                        break;
                    }

                    boolean wildcard = column.equals("*");
                    if (wildcard || column.equalsIgnoreCase(table.getColumnName(index))) {
                        matcher.reset(value.getStringValue(index));
                        if (matcher.find()) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) return false;
            }
            return true;  // If we get here we matched them all
        }

    }

    private void initComponents(Frame owner) {

        // All this to have tool tip text!
        table = new JTable() {
            @Override
            public String getToolTipText(MouseEvent e) {
                Point p = e.getPoint();
                int rowIndex = rowAtPoint(p);
                int colIndex = columnAtPoint(p);
                int realColumnIndex = convertColumnIndexToModel(colIndex);
                if (realColumnIndex > 0) {
                    return getModel().getValueAt(rowIndex, realColumnIndex).toString();
                }
                return null;
            }
        };

        //======== outer content pane ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        JPanel dialogPane = new JPanel();
        dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        dialogPane.setLayout(new BorderLayout());

        //======== main content panel ========
        JPanel contentPanel = new JPanel();
        contentPanel.setLayout(new BorderLayout(0, 10));

        //======== scrollPane1 ========
        JScrollPane scrollPane1 = new JScrollPane();
        scrollPane1.setViewportView(table);
        contentPanel.add(scrollPane1, BorderLayout.CENTER);

        //---- Filter  panel ----
        JPanel filterPanel = new JPanel();
        filterPanel.setLayout(new JideBoxLayout(filterPanel, JideBoxLayout.X_AXIS, 5));
        JLabel filterLabel = new JLabel("Filter:");
        final String filterToolTip = "Enter multiple filter strings separated by commas.  e.g.  GM12878, ChipSeq";
        filterLabel.setToolTipText(filterToolTip);
        filterPanel.add(filterLabel, JideBoxLayout.FIX);
        filterTextField = new JTextField();
        filterTextField.setToolTipText(filterToolTip);
        filterPanel.add(filterTextField, JideBoxLayout.VARY);

        rowCountLabel = new JLabel();
        rowCountLabel.setHorizontalAlignment(JLabel.RIGHT);
        JPanel rowCountPanel = new JPanel();
        rowCountPanel.setLayout(new JideBoxLayout(rowCountPanel, JideBoxLayout.X_AXIS, 0));
        rowCountPanel.setPreferredSize(new Dimension(100, 28));
        rowCountPanel.add(rowCountLabel, JideBoxLayout.VARY);
        filterPanel.add(rowCountPanel, JideBoxLayout.FIX);

        contentPanel.add(filterPanel, BorderLayout.NORTH);

        dialogPane.add(contentPanel, BorderLayout.CENTER);

        //======== buttonBar ========
        JPanel buttonBar = new JPanel();
        buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
        buttonBar.setLayout(new GridBagLayout());
        ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
        ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

        //---- okButton ----
        JButton okButton = new JButton("OK");
        getRootPane().setDefaultButton(okButton);
        okButton.addActionListener(e -> loadButtonActionPerformed(e));
        buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 5), 0, 0));

        //---- cancelButton ----
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> cancelButtonActionPerformed(e));

        buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));

        dialogPane.add(buttonBar, BorderLayout.SOUTH);

        contentPane.add(dialogPane, BorderLayout.CENTER);

        Rectangle ownerBounds = owner.getBounds();
        setSize(ownerBounds.width, 620);
        setLocationRelativeTo(getOwner());
    }


}
