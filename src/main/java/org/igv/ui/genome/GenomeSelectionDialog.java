/*
 * Created by JFormDesigner on Thu Oct 31 22:31:02 EDT 2013
 */

package org.igv.ui.genome;

import org.igv.DirectoryManager;
import org.igv.Globals;
import org.igv.event.GenomeResetEvent;
import org.igv.event.IGVEventBus;
import org.igv.feature.genome.GenomeDownloadUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.feature.genome.load.HubGenomeLoader;
import org.igv.feature.genome.load.JsonGenomeLoader;
import org.igv.feature.genome.load.TrackConfig;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ucsc.hub.Hub;
import org.igv.ucsc.hub.HubParser;
import org.igv.ui.util.MessageUtils;
import org.igv.util.LongRunningTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Jim Robinson
 */
public class GenomeSelectionDialog extends org.igv.ui.IGVDialog {

    private static Logger log = LogManager.getLogger(GenomeSelectionDialog.class);
    private static NumberFormatter numberFormatter = new NumberFormatter();

    private JTable table;
    private JTextField filterTextField;
    JPanel downloadSequencePanel;
    ButtonGroup downloadSequenceGroup;
    private JRadioButton downloadSequenceRB;
    private JRadioButton remoteSequenceRB;

    JPanel downloadAnnotationsPanel;
    ButtonGroup downloadAnnotationsGroup;
    private JRadioButton downloadAnnotationsRB;
    private JRadioButton remoteAnnotationsRB;
    private List<GenomeListItem> allListItems;
    private DefaultListModel genomeListModel;

    private GenomeTableModel model;
    private boolean canceled;

    public GenomeSelectionDialog(Frame owner, GenomeTableModel model) {
        super(owner);
        this.model = model;
        setModal(true);
        initComponents();
        init(model);
    }

    private void init(final GenomeTableModel model) {
        setModal(true);
        table.setRowSelectionAllowed(true);
        table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        table.setAutoCreateRowSorter(true);
        table.setModel(model);
        table.setRowSorter(model.getSorter());


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


        RowFilter<GenomeTableModel, Object> rf = null;
        //If current expression doesn't parse, don't update.
        try {
            rf = new RegexFilter(filterTextField.getText());
        } catch (java.util.regex.PatternSyntaxException e) {
            return;
        }
        model.getSorter().setRowFilter(rf);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
        IGVEventBus.getInstance().post(new GenomeResetEvent());
    }

    public boolean isCanceled() {
        return canceled;
    }

    private GenomeListItem getSelectedRecord() {
        int idx = table.getSelectedRow();
        if (idx < 0) {
            return null;
        } else {
            int modelIndex = table.convertRowIndexToModel(idx);
            return model.getRecords().get(modelIndex);
        }
    }

    private void initComponents() {

        //======== this ========
        setModal(true);
        setTitle("Genome Selector");

        JPanel outerPanel = new JPanel();
        outerPanel.setBorder(new EmptyBorder(20, 20, 20, 20));
        outerPanel.setLayout(new BorderLayout(0, 20));
        getContentPane().add(outerPanel);

        //======== contentPanel ========
        JPanel contentPanel = new JPanel();
        contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

        // Header message
        JLabel headerMessage = new JLabel();
        headerMessage.setText("Enter one or more search terms to filter genome list.");
        headerMessage.setRequestFocusEnabled(false);
        outerPanel.add(headerMessage, BorderLayout.NORTH);

        // Filter panel
        JPanel filterPanel = new JPanel();
        filterPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
        final String filterToolTip = "Enter multiple filter strings separated by spaces.";

        filterTextField = new JTextField();
        filterTextField.setMinimumSize(new Dimension(100, 25));
        filterTextField.setPreferredSize(new Dimension(900, 25));
        filterTextField.setToolTipText(filterToolTip);

        JLabel filterLabel = new JLabel("Filter:");
        filterLabel.setLabelFor(filterTextField);
        filterLabel.setRequestFocusEnabled(false);
        filterPanel.add(filterLabel);
        filterLabel.setToolTipText(filterToolTip);

        filterPanel.add(filterTextField);
        contentPanel.add(filterPanel);


        // Genome table
        table = new JTable();
        Font headerFont = table.getTableHeader().getFont();
        Font boldHeaderFont = headerFont.deriveFont(Font.BOLD);
        table.getTableHeader().setFont(boldHeaderFont);
        JScrollPane scrollPane1 = new JScrollPane();
        scrollPane1.setViewportView(table);
        contentPanel.add(scrollPane1);

        // Download buttons
        JPanel p1 = new JPanel();
        p1.setLayout(new GridLayout(2, 0));
        p1.setPreferredSize(new Dimension(500, 50));
        p1.setMinimumSize(new Dimension(500, 50));
        p1.setMaximumSize(new Dimension(1000, 50));

        downloadSequenceRB = new JRadioButton("Download sequence");
        remoteSequenceRB = new JRadioButton("Use remote sequence");
        remoteSequenceRB.setSelected(true);
        downloadSequenceGroup = new ButtonGroup();
        downloadSequenceGroup.add(downloadSequenceRB);
        downloadSequenceGroup.add(remoteSequenceRB);
        downloadSequencePanel = new JPanel();
        FlowLayout layout = new FlowLayout(FlowLayout.LEFT);
        downloadSequencePanel.setLayout(layout);
        downloadSequencePanel.add(downloadSequenceRB);
        downloadSequencePanel.add(remoteSequenceRB);
        p1.add(downloadSequencePanel);

        downloadAnnotationsRB = new JRadioButton("Download annotations");
        remoteAnnotationsRB = new JRadioButton("Use remote annotations");
        remoteAnnotationsRB.setSelected(true);
        downloadAnnotationsGroup = new ButtonGroup();
        downloadAnnotationsGroup.add(downloadAnnotationsRB);
        downloadAnnotationsGroup.add(remoteAnnotationsRB);
        downloadAnnotationsPanel = new JPanel();
        FlowLayout layout1 = new FlowLayout(FlowLayout.LEFT);
        downloadAnnotationsPanel.setLayout(layout1);
        downloadAnnotationsPanel.add(downloadAnnotationsRB);
        downloadAnnotationsPanel.add(remoteAnnotationsRB);
        p1.add(downloadAnnotationsPanel);
        contentPanel.add(p1);

        outerPanel.add(contentPanel, BorderLayout.CENTER);

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

        outerPanel.add(buttonBar, BorderLayout.SOUTH);

        if (this.model.getColumnCount() == 1) {
            setSize(500, 620);
        } else {
            setSize(1000, 620);
        }
        setLocationRelativeTo(getOwner());
    }

    private void loadButtonActionPerformed(ActionEvent evt) {
        canceled = false;
        try {
            GenomeListItem rec = getSelectedRecord();

            if (rec != null) {

                final String url = rec.getAttributeValue("url");
                final String id = rec.getAttributeValue("accession");

                Runnable showDialog = () -> {
                    try {

                        GenomeConfig config;
                        if (url != null && url.endsWith(".json")) {
                            config = (new JsonGenomeLoader(url)).loadGenomeConfig();
                        } else {
                            String accession = rec.getAttributeValue("accession");
                            String hubURL = HubGenomeLoader.convertToHubURL(accession);
                            Hub hub = HubParser.loadHub(hubURL);
                            config = hub.getGenomeConfigs().get(0);
                            config.setHubs(Arrays.asList(hubURL));
                        }

                        config.setName(rec.getAttributeValue("common name"));

                        // If config has a hub,  allow changing default annotation.
                        if (config.getHubs() != null && config.getHubs().size() > 0) {

                            List<TrackConfig> selectedTracks = GenomeManager.selectAnnotationTracks(config, GenomeManager.SELECT_ANNOTATIONS_MESSAGE);
                            if (selectedTracks != null && selectedTracks.size() > 0) {
                                config.setTracks(selectedTracks);
                            }
                        }

                        File localFile = GenomeDownloadUtils.downloadGenome(config,
                                downloadSequenceRB.isSelected(),
                                downloadAnnotationsRB.isSelected());

                        if (localFile != null) {
                            GenomeManager.getInstance().loadGenome(localFile.getAbsolutePath());
                        } else {
                            GenomeManager.getInstance().loadGenome(url);
                        }

                        // Legacy cleanup
                        removeDotGenomeFile(id);

                    } catch (IOException e) {
                        MessageUtils.showErrorMessage("Error loading genome " + url, e);
                        log.error("Error loading genome " + url, e);
                    }
                };

                if (SwingUtilities.isEventDispatchThread()) {
                    LongRunningTask.submit(showDialog);
                } else {
                    showDialog.run();
                }

            }
        } catch (Exception ex) {
            log.error(ex);
            MessageUtils.showMessage("Error loading genome: " + ex.getMessage());
        }

        setVisible(false);
    }

    private static void removeDotGenomeFile(String id) {
        try {
            File dotGenomeFile = new File(DirectoryManager.getGenomeCacheDirectory(), id + ".genome");
            if (dotGenomeFile.exists()) {
                dotGenomeFile.delete();
            }
        } catch (Exception e) {
            // If anything goes wrong, just log it, this cleanup is not essential
            log.error("Error deleting .genome file", e);
        }
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
