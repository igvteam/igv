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

package org.broad.igv.ui.table;

import com.jidesoft.swing.JideBoxLayout;
import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.ParsingUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Jim Robinson
 */
public class SearchableTableDialog extends org.broad.igv.ui.IGVDialog {

    private static Logger log = LogManager.getLogger(SearchableTableDialog.class);
    private static NumberFormatter numberFormatter = new NumberFormatter();

    private JTable table;
    private JTextField filterTextField;
    private JLabel rowCountLabel;
    SearchableTableModel model;
    private boolean canceled;

    public SearchableTableDialog(Frame owner, SearchableTableModel model) {
        super(owner);
        this.model = model;
        setModal(true);
        initComponents();
        init(model);
    }

    private void init(final SearchableTableModel model) {
        setModal(true);
        //setTitle("Encode Production Data");

        table.setRowSelectionAllowed(true);
        table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        table.setAutoCreateRowSorter(true);
        table.setModel(model);
        table.setRowSorter(model.getSorter());

        try {
            rowCountLabel.setText(numberFormatter.valueToString(table.getRowCount()) + " rows");
        } catch (ParseException e) {
            rowCountLabel.setText("" + table.getRowCount() + " rows");
        }

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

    /**
     * Update the row filter regular expression from the expression in
     * the text box.
     */
    private void updateFilter() {


        RowFilter<SearchableTableModel, Object> rf = null;
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
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
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

    public SearchableTableRecord getSelectedRecord() {
        int idx = table.getSelectedRow();
        if(idx < 0) {
            return  null;
        } else {
            int modelIndex = table.convertRowIndexToModel(idx);
            return model.getRecords().get(modelIndex);
        }
    }

    public List<SearchableTableRecord> getSelectedRecords() {
        return Arrays.stream(table.getSelectedRows())
                .mapToObj(table::convertRowIndexToModel)
                .map(model.getRecords()::get)
                .collect(Collectors.toList());
    }

    private void initComponents() {


        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== contentPanel ========
        JPanel contentPanel = new JPanel();
        contentPanel.setLayout(new BorderLayout(0, 10));

        table = new JTable();
        filterTextField = new JTextField();
        rowCountLabel = new JLabel();
        rowCountLabel.setHorizontalAlignment(JLabel.RIGHT);

        Font headerFont = table.getTableHeader().getFont();
        Font boldHeaderFont = headerFont.deriveFont(Font.BOLD);
        table.getTableHeader().setFont(boldHeaderFont);

        //======== dialogPane ========
        JPanel dialogPane = new JPanel();
        dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        dialogPane.setLayout(new BorderLayout());

        //======== scrollPane1 ========
        JScrollPane scrollPane1 = new JScrollPane();
        scrollPane1.setViewportView(table);
        contentPanel.add(scrollPane1, BorderLayout.CENTER);

        //======== filterPanel ========
        JPanel filterPanel = new JPanel();
        filterPanel.setLayout(new JideBoxLayout(filterPanel, JideBoxLayout.X_AXIS, 5));
        filterPanel.add(filterTextField, JideBoxLayout.VARY);

        //---- filterLabel ----
        JLabel filterLabel = new JLabel();
        filterLabel.setText("Filter:");
        filterPanel.add(filterLabel, JideBoxLayout.FIX);
        final String filterToolTip = "Enter multiple filter strings separated by commas.  e.g.  GM12878, ChipSeq";
        filterLabel.setToolTipText(filterToolTip);
        filterTextField.setToolTipText(filterToolTip);


        JPanel sillyPanel = new JPanel();
        sillyPanel.setLayout(new JideBoxLayout(sillyPanel, JideBoxLayout.X_AXIS, 0));
        sillyPanel.setPreferredSize(new Dimension(100, 28));
        sillyPanel.add(rowCountLabel, JideBoxLayout.VARY);
        filterPanel.add(sillyPanel, JideBoxLayout.FIX);

        contentPanel.add(filterPanel, BorderLayout.NORTH);
        dialogPane.add(contentPanel, BorderLayout.CENTER);

        //======== buttonBar ========

        JPanel buttonBar = new JPanel();
        buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
        buttonBar.setLayout(new GridBagLayout());
        ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
        ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

        //---- okButton ----
        JButton okButton = new JButton();
        getRootPane().setDefaultButton(okButton);
        okButton.setText("Load");
        okButton.addActionListener(e -> loadButtonActionPerformed(e));
        buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 5), 0, 0));

        //---- cancelButton ----
        JButton cancelButton = new JButton();
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(e -> cancelButtonActionPerformed(e));
        buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));

        dialogPane.add(buttonBar, BorderLayout.SOUTH);

        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(1000, 620);
        setLocationRelativeTo(getOwner());
    }

    private class RegexFilter extends RowFilter {

        List<Matcher> matchers;


        RegexFilter(String text) {
            if (text == null) {
                throw new IllegalArgumentException("Pattern must be non-null");
            }

            matchers = Arrays.stream(Globals.whitespacePattern.split(text))
                    .map(t -> {
                        String value = t.trim();
                        return Pattern.compile("(?i)" + value).matcher("");
                    })
                    .collect(Collectors.toList());
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
            return matchers.stream()
                    .allMatch(entry -> {
                        Matcher matcher = entry;

                        return IntStream.range(0, table.getColumnCount())
                                .anyMatch(index -> {
                                    matcher.reset(table.getColumnName(index).toLowerCase());
                                    if (matcher.find()) {
                                        return true;
                                    }

                                    return matcher.reset(value.getStringValue(index)).find();
                                });
                    });
        }

    }

}
