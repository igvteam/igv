package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.genome.GenomeTableModel;

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

public class HubSelectionDialog extends JDialog {

    private static final long serialVersionUID = 1L;

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

        // Filter track hubs to those supporting the superset of IGV loaded + UCSC DB genomes
        String genomeID = GenomeManager.getInstance().getCurrentGenome().getUCSCId();

        List<HubDescriptor> allHubDescriptors = HubRegistry.getAllHubsForGenome(genomeID);

        final List<HubDescriptor> allSelectedHubs = HubRegistry.getSelectedHubsForGenome(genomeID);

        this.hubDescriptors = allHubDescriptors.stream()
                .filter(hd -> hd.getDbList() != null && hd.getDbList().contains(genomeID))
                .collect(Collectors.toList());

        if (this.hubDescriptors.isEmpty()) {
            JOptionPane.showMessageDialog(this, "No hubs available");
            return;
        }

        initComponents(allSelectedHubs);
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

        JLabel titleLabel = new JLabel("Track hubs for genome: " + GenomeManager.getInstance().getCurrentGenome().getDisplayName());
        titleLabel.setFont(titleLabel.getFont().deriveFont(Font.BOLD, 16));
        headerPanel.add(titleLabel);
        headerPanel.add(Box.createVerticalStrut(20));

        JLabel headerMessage = new JLabel();
        headerMessage.setText("<html>Hub list from <a href='https://genome.ucsc.edu/cgi-bin/hgHubConnect'>https://genome.ucsc.edu/cgi-bin/hgHubConnect</a></html>");
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

