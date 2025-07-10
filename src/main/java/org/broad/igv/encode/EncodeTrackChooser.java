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
import java.util.stream.Collectors;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.NumberFormatter;

import com.jidesoft.swing.JideBoxLayout;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.BrowseEncodeAction;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.Pair;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

/**
 * @author Jim Robinson
 */
public class EncodeTrackChooser extends org.broad.igv.ui.IGVDialog {

    private static Logger log = LogManager.getLogger(EncodeTrackChooser.class);

    private static Map<String, EncodeTrackChooser> instanceMap = Collections.synchronizedMap(new HashMap<>());
    private static NumberFormatter numberFormatter = new NumberFormatter();

    private static String ENCODE_HOST = "https://www.encodeproject.org";
    private static Set<String> filteredColumns = new HashSet(Arrays.asList("ID", "Assembly", "HREF", "path"));

    private static List<String> filteredExtensions = Arrays.asList("tsv", "tsv.gz");

    private static Map<String, String> speciesNames = Map.of(
            "ce10", "Caenorhabditis elegans",
            "ce11", "Caenorhabditis elegans",
            "dm3", "Drosophila melanogaster",
            "dm6", "Drosophila melanogaster",
            "GRCh38", "Homo sapiens",
            "hg19", "Homo sapiens",
            "mm10", "Mus musculus",
            "mm9", "Mus musculus"
    );

    static HashSet<String> ucscSupportedGenomes = new HashSet<>(Arrays.asList("hg19", "mm9"));
    static HashSet<String> supportedGenomes = new HashSet<>(
            Arrays.asList("ce10", "ce11", "dm3", "dm6", "GRCh38", "hg19", "mm10", "mm9"));



    JTable table;
    JTextField filterTextField;
    JLabel rowCountLabel;
    EncodeTableModel model;
    private boolean canceled;


    /**
     * Return a new or cached instance of a track chooser for the given genome and type.
     *
     * @param genomeId
     * @param type
     * @return
     * @throws IOException
     */
    public synchronized static EncodeTrackChooser getInstance(String genomeId, BrowseEncodeAction.Type type) throws IOException {

        String encodeGenomeId = getEncodeGenomeID(genomeId);
        String key = encodeGenomeId + type.toString();
        EncodeTrackChooser instance = instanceMap.get(key);
        if (instance == null) {
            Pair<List<String>, List<EncodeFileRecord>> records = getEncodeFileRecords(encodeGenomeId, type);
            if (records == null) {
                return null;
            }
            Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
            final List<String> headings = records.getFirst();
            final List<EncodeFileRecord> rows = records.getSecond();
            final String title = getDialogTitle(genomeId, type);
            instance = new EncodeTrackChooser(parent, new EncodeTableModel(headings, rows), title);
            instanceMap.put(key, instance);
        }

        return instance;
    }

    private static String getDialogTitle(String genomeId, BrowseEncodeAction.Type type) {

        if (type == BrowseEncodeAction.Type.UCSC) {
            return "ENCODE data hosted at UCSC (2012)";
        } else {
            switch (type) {
                case SIGNALS_CHIP:
                    return "ENCODE CHiP Seq - Signals";
                case SIGNALS_OTHER:
                    return "ENCODE Non CHiP Data - Signals";
                default:
                    return "ENCODE";
            }
        }
    }

    public static boolean genomeSupportedUCSC(String genomeId) {
        return genomeId != null && ucscSupportedGenomes.contains(getEncodeGenomeID(genomeId));
    }

    public static boolean genomeSupported(String genomeId) {
        return genomeId != null && supportedGenomes.contains(getEncodeGenomeID(genomeId));
    }


    private static String getEncodeGenomeID(String genomeId) {
        switch (genomeId) {
            case "hg38":
            case "hg38_1kg":
                return "GRCh38";
            case "b37":
            case "1kg_v37":
                return "hg19";
            default:
                return genomeId;
        }

    }

    private static Pair<List<String>, List<EncodeFileRecord>> getEncodeFileRecords(String genomeId, BrowseEncodeAction.Type type) throws IOException {

        try (InputStream is = getStreamFor(genomeId, type)) {
            if (is == null) {
                return null;
            }
            return parseRecords(is, type, genomeId);
        }
    }

    private static InputStream getStreamFor(String genomeId, BrowseEncodeAction.Type type) throws IOException {
        if (type == BrowseEncodeAction.Type.UCSC) {
            return EncodeTrackChooser.class.getResourceAsStream("encode." + genomeId + ".txt");
        } else {
            String root = PreferencesManager.getPreferences().get(Constants.ENCODE_FILELIST_URL) + genomeId + ".";
            String url = null;
            switch (type) {
                case SIGNALS_CHIP:
                    url = root + "signals.chip.txt.gz";
                    break;
                case SIGNALS_OTHER:
                    url = root + "signals.other.txt.gz";
                    break;
                case OTHER:
                    url = root + "other.txt.gz";
                    break;
            }
            if (url == null) {
                throw new RuntimeException("Unknown encode data collection type: " + type);
            }
            return ParsingUtils.openInputStream(url);
        }
    }

    private static Pair parseRecords(InputStream is, BrowseEncodeAction.Type type, String genomeId) throws IOException {

        BufferedReader reader = new BufferedReader(new InputStreamReader(is));

        String[] headers = Globals.tabPattern.split(reader.readLine());

        int pathColumn = type == BrowseEncodeAction.Type.UCSC ? 0 : Arrays.asList(headers).indexOf("HREF");

        List<EncodeFileRecord> records = new ArrayList<>(20000);
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            if (!nextLine.startsWith("#")) {

                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                String path = type == BrowseEncodeAction.Type.UCSC ? tokens[pathColumn] : ENCODE_HOST + tokens[pathColumn];

                if(filteredExtensions.stream().anyMatch(e -> path.endsWith(e))) {
                    continue;
                }

                Map<String, String> attributes = new LinkedHashMap<>();
                for (int i = 0; i < headers.length; i++) {
                    String value = i < tokens.length ? tokens[i] : "";
                    if (value.length() > 0) {
                        attributes.put(headers[i], shortenField(value, genomeId));
                    }
                }
                final EncodeFileRecord record = new EncodeFileRecord(path, attributes);
                records.add(record);

            }
        }

        List<String> filteredHeaders = Arrays.stream(headers).filter(h -> !filteredColumns.contains(h)).collect(Collectors.toList());

        return new Pair(filteredHeaders, records);
    }

    private static String shortenField(String value, String genomeId) {
        String species = speciesNames.get(genomeId);
        return species == null ?
                value :
                value.replace("(" + species + ")", "").replace(species, "").trim();
    }


    private EncodeTrackChooser(Frame owner, EncodeTableModel model, String title) {
        super(owner);
        setTitle(title);
        this.model = model;
        setModal(true);
        initComponents(owner);
        init(model);
    }

    private void init(final EncodeTableModel model) {
        setModal(true);

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
                record.setSelected(false);   // Prevent loading twice
            }
        }

        return selectedRecords;
    }

    private class RegexFilter extends RowFilter {

        List<Pair<String, Matcher>> matchers;

        RegexFilter(String text) {

            if (text == null) {
                throw new IllegalArgumentException("Pattern must be non-null");
            }
            matchers = new ArrayList<Pair<String, Matcher>>();
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
                java.awt.Point p = e.getPoint();
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
        JButton okButton = new JButton("Load");
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


    /**
     * Main function for testing only
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        getInstance("hg19", BrowseEncodeAction.Type.UCSC).setVisible(true);
    }

}
