/*
 * Created by JFormDesigner on Thu Oct 31 22:31:02 EDT 2013
 */

package org.broad.igv.util.encode;

import java.awt.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.NumberFormatter;

import com.jidesoft.swing.JideBoxLayout;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ResourceLocator;

/**
 * @author Jim Robinson
 */
public class EncodeFileBrowser extends JDialog {

    private static Logger log = Logger.getLogger(EncodeFileBrowser.class);

    private static EncodeFileBrowser theInstance;
    private static NumberFormatter numberFormatter = new NumberFormatter();

    private JButton cancelButton;
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JTable table;
    private JPanel filterPanel;
    private JLabel filterLabel;
    private JTextField filterTextField;
    private JLabel rowCountLabel;
    private JPanel buttonBar;
    private JButton okButton;

    EncodeTableModel model;
    private boolean canceled;


    public synchronized static EncodeFileBrowser getInstance() throws IOException {

        if (theInstance == null) {
            List<EncodeFileRecord> records = getEncodeFileRecords();
            Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
            theInstance = new EncodeFileBrowser(parent, new EncodeTableModel(records));
        }

        return theInstance;
    }

    private static List<EncodeFileRecord> getEncodeFileRecords() throws IOException {

        InputStream is = null;

        try {

            is = EncodeFileBrowser.class.getResourceAsStream("encode.txt");
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String[] headers = Globals.tabPattern.split(reader.readLine());
            List<EncodeFileRecord> records = new ArrayList<EncodeFileRecord>(20000);
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (!nextLine.startsWith("#")) {

                    String[] tokens = Globals.tabPattern.split(nextLine, -1);
                    String path = tokens[0];

                    Map<String, String> attributes = new HashMap<String, String>();
                    for (int i = 0; i < headers.length; i++) {
                        String value = tokens[i];
                        if (value.length() > 0) {
                            attributes.put(headers[i], value);
                        }
                    }
                    records.add(new EncodeFileRecord(path, attributes));

                }

            }
            return records;
        } finally {
            if (is != null) is.close();
        }
    }


    private EncodeFileBrowser(Frame owner, EncodeTableModel model) {
        super(owner);
        this.model = model;
        setModal(true);
        initComponents();
        init(model);
    }

    private void init(final EncodeTableModel model) {
        setModal(true);
        setTitle("Encode Production Data");

        table.setAutoCreateRowSorter(true);
        table.setModel(model);
        table.setRowSorter(model.getSorter());
        try {
            rowCountLabel.setText(numberFormatter.valueToString(table.getRowCount()) + " rows");
        } catch (ParseException e) {

        }

        table.setRowSelectionAllowed(false);
        table.setColumnSelectionAllowed(false);

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


        RowFilter<EncodeTableModel, Object> rf = null;
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

    private Set<String> getLoadedPaths() {

        if (!IGV.hasInstance()) return new HashSet<String>();

        Collection<ResourceLocator> locators = IGV.getInstance().getDataResourceLocators();
        HashSet<String> loadedPaths = new HashSet<String>(locators.size());
        for (ResourceLocator locator : locators) {
            loadedPaths.add(locator.getPath());
        }
        return loadedPaths;
    }

    /**
     * @return the list of VISIBLE selected records.  Filtered records are not returned even if record.selected == true
     * @throws IOException
     */
    public List<EncodeFileRecord> getSelectedRecords() throws IOException {

        List<EncodeFileRecord> selectedRecords = new ArrayList<EncodeFileRecord>();
        List<EncodeFileRecord> allRecords = model.getRecords();

        int rowCount = table.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            int modelIdx = table.convertRowIndexToModel(i);
            EncodeFileRecord record = allRecords.get(modelIdx);
            if (record.isSelected()) {
                selectedRecords.add(record);
            }
        }

        return selectedRecords;
    }

    private class RegexFilter extends RowFilter {

        Map<String, Matcher> matchers;

        RegexFilter(String text) {

            if (text == null) {
                throw new IllegalArgumentException("Pattern must be non-null");
            }
            matchers = new HashMap<String, Matcher>();
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
                    }
                    else {
                        value = kv[0];  // Value is column name until more input is entered
                    }
                }

                matchers.put(column, Pattern.compile("(?i)" + value).matcher(""));
            }

        }


        @Override
        public boolean include(Entry value) {

            for (Map.Entry<String, Matcher> entry : matchers.entrySet()) {
                String column = entry.getKey();
                Matcher matcher = entry.getValue();

                int nColumns = table.getColumnCount();
                // First column is checkbox

                boolean found = false;
                for (int index = 1; index < nColumns; index++) {

                    // Include column headings in search
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
                if (!found) return false; // End of column loop.  Must find a match for all matchers
            }
            return true;  // If we get here we matched them all
        }
    }


    private void initComponents() {

        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        table = new JTable();
        filterPanel = new JPanel();
        filterLabel = new JLabel();
        filterTextField = new JTextField();
        rowCountLabel = new JLabel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        getRootPane().setDefaultButton(okButton);

        final String filterToolTip = "Enter multiple filter strings separated by commas.  e.g.  GM12878, ChipSeq";
        filterLabel.setToolTipText(filterToolTip);
        filterTextField.setToolTipText(filterToolTip);

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========

        dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        dialogPane.setLayout(new BorderLayout());

        //======== contentPanel ========

        contentPanel.setLayout(new BorderLayout(0, 10));

        //======== scrollPane1 ========

        scrollPane1.setViewportView(table);

        contentPanel.add(scrollPane1, BorderLayout.CENTER);

        //======== panel1 ========

        filterPanel.setLayout(new JideBoxLayout(filterPanel, JideBoxLayout.X_AXIS, 5));

        //---- label1 ----
        filterLabel.setText("Filter:");
        filterPanel.add(filterLabel, JideBoxLayout.FIX);

        //---- filterTextField ----
        filterPanel.add(filterTextField, JideBoxLayout.VARY);

        rowCountLabel.setHorizontalAlignment(JLabel.RIGHT);
        JPanel sillyPanel = new JPanel();
        sillyPanel.setLayout(new JideBoxLayout(sillyPanel, JideBoxLayout.X_AXIS, 0));
        sillyPanel.setPreferredSize(new Dimension(100, 28));
        sillyPanel.add(rowCountLabel, JideBoxLayout.VARY);

        filterPanel.add(sillyPanel, JideBoxLayout.FIX);

        contentPanel.add(filterPanel, BorderLayout.NORTH);

        dialogPane.add(contentPanel, BorderLayout.CENTER);

        //======== buttonBar ========

        buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
        buttonBar.setLayout(new GridBagLayout());
        ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
        ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

        //---- okButton ----
        okButton.setText("Load");
        okButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                loadButtonActionPerformed(e);
            }
        });
        buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 5), 0, 0));

        //---- cancelButton ----
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                cancelButtonActionPerformed(e);
            }
        });
        buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));

        dialogPane.add(buttonBar, BorderLayout.SOUTH);

        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(700, 620);
        setLocationRelativeTo(getOwner());
    }


    public static void main(String[] args) throws IOException {
        getInstance().setVisible(true);
    }

}
