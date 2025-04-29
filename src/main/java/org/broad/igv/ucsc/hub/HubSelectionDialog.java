package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.genome.GenomeTableModel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class HubSelectionDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    private static Map<String, HubSelectionDialog> instances = new java.util.HashMap<>();

    List<HubDescriptor> hubDescriptors;
    private boolean canceled;
    private JTable table;
    private JTextField filterTextField;

    public HubSelectionDialog(Frame owner) {

        super(owner, "Select Track Hubs", true);

        setSize(owner.getWidth() - 20, 600);
        if (owner != null) {
            setLocationRelativeTo(owner);
        }


        List<HubDescriptor> allHubDescriptors = HubRegistry.getAllHubs();
        final List<HubDescriptor> allSelectedHubs = HubRegistry.getAllSelectedHubs();

        // Filter track hubs to those supporting the superset of IGV loaded + UCSC DB genomes
        allHubDescriptors = filterHubs(allHubDescriptors, allSelectedHubs);

        if (allHubDescriptors != null && !allHubDescriptors.isEmpty()) {
            this.hubDescriptors = allHubDescriptors;
        } else {
            this.hubDescriptors = new ArrayList<>();
        }

        if (this.hubDescriptors.isEmpty()) {
            JOptionPane.showMessageDialog(this, "No hubs available");
            return;
        }

        initComponents(allSelectedHubs);

        // Preset filter to currently loaded genome
        String genomeID = GenomeManager.getInstance().getCurrentGenome().getId();
        filterTextField.setText(genomeID);

    }

    private List<HubDescriptor> filterHubs(List<HubDescriptor> hubDescriptors, List<HubDescriptor> allSelectedHubs) {

        Set<String> allSelectedURLs = allSelectedHubs.stream().map(HubDescriptor::getUrl).collect(Collectors.toSet());

        return hubDescriptors.stream().filter(h -> {
            Set<String> genomeIDs = getLoadedGenomeIDs();
            genomeIDs.addAll(HubRegistry.getUcscGenomeIDs());
            String[] hubGenomes = h.getDbList().split(",");
            if (hubGenomes == null || hubGenomes.length == 0) {
                return true;
            }
            if (allSelectedURLs.contains(h.getUrl())) {
                return true;
            }
            for (String genome : hubGenomes) {
                if (genomeIDs.contains(genome)) {
                    return true;
                }
            }
            return false;
        }).toList();
    }

    /**
     * Get a list of genome IDs from the genome list manager.
     * @return
     */
    private Set<String> getLoadedGenomeIDs() {
        try {
            Collection<GenomeListItem> genomeListItems = GenomeListManager.getInstance().getGenomeItemMap().values();
            HashSet<String> ids = new HashSet<>(genomeListItems.stream().map(GenomeListItem::getId).collect(Collectors.toList()));
            ids.add("hg38");
            ids.add("hg19");
            return ids;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    public List<HubDescriptor> getSelectedHubs() {
        return this.hubDescriptors.stream().filter(HubDescriptor::isSelected).collect(Collectors.toList());
    }


    private void initComponents(Collection<HubDescriptor> loadedHubs) {

        JPanel outerPanel = new JPanel();
        outerPanel.setBorder(new EmptyBorder(10, 10, 10, 10));
        outerPanel.setLayout(new BorderLayout(0, 20));
        add(outerPanel);

        // Header message
        JPanel headerPanel = new JPanel();
        headerPanel.setLayout(new BoxLayout(headerPanel, BoxLayout.Y_AXIS));
        JLabel headerMessage = new JLabel("Public track hubs from https://genome.ucsc.edu/goldenpath/help/api.html");
        headerPanel.add(headerMessage);
        outerPanel.add(headerPanel, BorderLayout.NORTH);

        //======== contentPanel ========
        JPanel contentPanel = new JPanel();
        contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
        outerPanel.add(contentPanel, BorderLayout.CENTER);

        // Filter panel
        JPanel filterPanel = new JPanel();
        filterPanel.setLayout(new BoxLayout(filterPanel, BoxLayout.X_AXIS));
        final String filterToolTip = "Enter multiple filter strings separated by spaces.";

        JLabel filterLabel = new JLabel("Filter:");
        filterLabel.setLabelFor(filterTextField);
        filterLabel.setRequestFocusEnabled(false);
        filterLabel.setToolTipText(filterToolTip);
        filterPanel.add(filterLabel);

        filterTextField = new JTextField();
        filterTextField.setToolTipText(filterToolTip);
        filterTextField.setMaximumSize(new Dimension(Integer.MAX_VALUE, filterTextField.getPreferredSize().height));
        filterTextField.setToolTipText("Type keywords to filter the track hubs, e.g. hg38. Use spaces to separate multiple keywords.");

        filterPanel.add(Box.createHorizontalStrut(5)); // Add spacing between label and text field
        filterPanel.add(filterTextField);

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

        contentPanel.add(filterPanel);


        // Table
        Set<String> loadedURLs = loadedHubs.stream().map(HubDescriptor::getUrl).collect(Collectors.toSet());
        for (int i = 0; i < this.hubDescriptors.size(); i++) {
            HubDescriptor hub = this.hubDescriptors.get(i);
            hub.setSelected(loadedURLs.contains(hub.getUrl()));
        }
        HubTableModel model = new HubTableModel(this.hubDescriptors);
        table = new HubTable(model);

        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
        contentPanel.add(scrollPane);


        // Button panel
        JPanel buttonPanel = new JPanel();
        ((FlowLayout) buttonPanel.getLayout()).setAlignment(FlowLayout.RIGHT);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> {
            canceled = true;
            setVisible(false);
            dispose();
        });

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> {
            canceled = false;
            setVisible(false);
            dispose();
        });

//        if (Globals.IS_MAC) {
//            buttonPanel.add(cancelButton);
//            buttonPanel.add(okButton);
//        } else {
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        //}
        outerPanel.add(buttonPanel, BorderLayout.SOUTH);
    }

    public boolean isCanceled() {
        return canceled;
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
        ((DefaultRowSorter) table.getRowSorter()).setRowFilter(rf);
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

