package org.igv.ucsc.hub;

import org.igv.Globals;
import org.igv.feature.genome.GenomeManager;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Single-select chooser over the UCSC public hub list.  Hubs already added, from this list or by URL, are shown
 * disabled.  The chosen hub is loaded by the caller exactly as a hub entered by URL.
 */
public class HubSelectionDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    List<HubDescriptor> hubDescriptors;
    private boolean canceled = true;
    private HubTable table;
    private JTextField filterTextField;
    private JCheckBox currentGenomeCheckBox;
    private Set<String> currentGenomeUrls;

    public HubSelectionDialog(Frame owner) {

        super(owner, "Select Public Hub", true);

        setSize(owner.getWidth() - 20, 600);
        if (owner != null) {
            setLocationRelativeTo(owner);
        }

        // UCSC public hubs for all genomes, with hubs for the current genome listed first.  Hubs already added
        // (from this list or by URL) are shown disabled.
        String genomeID = GenomeManager.getInstance().getCurrentGenome().getUCSCId();

        currentGenomeUrls = HubRegistry.getAllHubsForGenome(genomeID).stream()
                .map(HubDescriptor::getUrl)
                .collect(Collectors.toSet());

        Set<String> registeredUrls = HubRegistry.getSelectedHubs().stream()
                .map(HubDescriptor::getUrl)
                .collect(Collectors.toSet());

        this.hubDescriptors = HubRegistry.getUcscPublicHubDescriptors().stream()
                .sorted(Comparator.comparing(hd -> currentGenomeUrls.contains(hd.getUrl()) ? 0 : 1))
                .collect(Collectors.toList());

        initComponents(registeredUrls);
    }


    /**
     * @return the chosen hub, or null if the dialog was canceled or no hub was selected.
     */
    public HubDescriptor getSelectedHub() {
        if (canceled) {
            return null;
        }
        int row = table.getSelectedRow();
        return row < 0 ? null : this.hubDescriptors.get(table.convertRowIndexToModel(row));
    }


    private void initComponents(Set<String> registeredUrls) {

        JPanel outerPanel = new JPanel();
        outerPanel.setBorder(new EmptyBorder(10, 10, 10, 10));
        outerPanel.setLayout(new BorderLayout(0, 20));
        add(outerPanel);

        // Header message
        JPanel headerPanel = new JPanel();
        headerPanel.setLayout(new BoxLayout(headerPanel, BoxLayout.Y_AXIS));

        JLabel titleLabel = new JLabel("UCSC public hubs");
        titleLabel.setFont(titleLabel.getFont().deriveFont(Font.BOLD, 16));
        headerPanel.add(titleLabel);
        headerPanel.add(Box.createVerticalStrut(20));

        JLabel headerMessage = new JLabel();
        headerMessage.setText("<html>Select a hub from the <a href='https://genome.ucsc.edu/cgi-bin/hgHubConnect'>" +
                "UCSC Track Hub Registry</a> to add to the 'Tracks' menu.  Hubs already added are greyed out.<br>" +
                "Hubs for the current genome (" + GenomeManager.getInstance().getCurrentGenome().getDisplayName() +
                ") are listed first; uncheck 'Current genome only' to see hubs for other genomes.</html>");
        headerMessage.setCursor(new Cursor(Cursor.HAND_CURSOR));
        headerMessage.addMouseListener(new java.awt.event.MouseAdapter() {
            @Override
            public void mouseClicked(java.awt.event.MouseEvent e) {
                try {
                    Desktop.getDesktop().browse(new java.net.URI("https://genome.ucsc.edu/cgi-bin/hgHubConnect"));
                } catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog(HubSelectionDialog.this,
                            "Unable to open the web browser. Please visit:\nhttps://genome.ucsc.edu/cgi-bin/hgHubConnect",
                            "Error Opening Browser",
                            JOptionPane.ERROR_MESSAGE);
                }
            }
        });
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

        currentGenomeCheckBox = new JCheckBox("Current genome only", true);
        currentGenomeCheckBox.setToolTipText("Show only hubs with tracks for the current genome");
        currentGenomeCheckBox.addActionListener(e -> updateFilter());
        filterPanel.add(Box.createHorizontalStrut(10));
        filterPanel.add(currentGenomeCheckBox);

        contentPanel.add(filterPanel);


        // Table
        HubTableModel model = new HubTableModel(this.hubDescriptors);
        table = new HubTable(model, registeredUrls);

        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
        contentPanel.add(scrollPane);

        updateFilter();


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

    /**
     * Update the row filter from the expression in the text box and the "current genome only" checkbox.
     */
    private void updateFilter() {

        RowFilter<HubTableModel, Integer> rf = null;
        //If current expression doesn't parse, don't update.
        try {
            rf = new RegexFilter(filterTextField.getText());
        } catch (java.util.regex.PatternSyntaxException e) {
            return;
        }
        if (currentGenomeCheckBox.isSelected()) {
            // The genome filter is the cheap one and rejects most rows, so it goes first
            rf = RowFilter.andFilter(Arrays.asList(new GenomeFilter(), rf));
        }
        ((DefaultRowSorter) table.getRowSorter()).setRowFilter(rf);
    }

    /**
     * Include only hubs listing the current genome.
     */
    private class GenomeFilter extends RowFilter<HubTableModel, Integer> {
        @Override
        public boolean include(Entry<? extends HubTableModel, ? extends Integer> entry) {
            return currentGenomeUrls.contains(entry.getModel().getRecords().get(entry.getIdentifier()).getUrl());
        }
    }

    private class RegexFilter extends RowFilter<HubTableModel, Integer> {

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
        public boolean include(Entry<? extends HubTableModel, ? extends Integer> value) {
            return matchers.stream()
                    .allMatch(entry -> {
                        Matcher matcher = entry;

                        // Model column 0 is the checkbox, not hub text
                        return IntStream.range(1, value.getValueCount())
                                .anyMatch(index -> {
                                    matcher.reset(value.getModel().getColumnName(index).toLowerCase());
                                    if (matcher.find()) {
                                        return true;
                                    }

                                    return matcher.reset(value.getStringValue(index)).find();
                                });
                    });
        }

    }
}
