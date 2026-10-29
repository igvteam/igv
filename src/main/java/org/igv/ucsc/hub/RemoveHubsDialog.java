package org.igv.ucsc.hub;

import org.igv.feature.genome.GenomeManager;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Checkbox list of the hubs in the registry (selected from the UCSC public list or added by URL), with a "Remove"
 * button.  Modeled on "Genomes > Remove Genomes...".
 */
public class RemoveHubsDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    private List<HubDescriptor> hubDescriptors;
    private boolean canceled = true;

    public RemoveHubsDialog(Frame owner) {

        super(owner, "Remove Hubs", true);

        setSize(owner.getWidth() - 20, 600);
        if (owner != null) {
            setLocationRelativeTo(owner);
        }

        // All registry hubs, for all genomes, with hubs for the current genome listed first.  The "selected" flag is
        // the checkbox state: checked hubs are removed.
        String genomeID = GenomeManager.getInstance().getCurrentGenome().getUCSCId();

        Set<String> currentGenomeUrls = HubRegistry.getSelectedHubsForGenome(genomeID).stream()
                .map(HubDescriptor::getUrl)
                .collect(Collectors.toSet());

        this.hubDescriptors = HubRegistry.getSelectedHubs().stream()
                .sorted(Comparator.comparing(hd -> currentGenomeUrls.contains(hd.getUrl()) ? 0 : 1))
                .collect(Collectors.toList());

        for (HubDescriptor hub : this.hubDescriptors) {
            hub.setSelected(false);
        }

        initComponents();
    }

    /**
     * @return the hubs checked for removal; empty if the dialog was canceled.
     */
    public List<HubDescriptor> getHubsToRemove() {
        if (canceled) {
            return Collections.emptyList();
        }
        return this.hubDescriptors.stream().filter(HubDescriptor::isSelected).collect(Collectors.toList());
    }

    private void initComponents() {

        JPanel outerPanel = new JPanel();
        outerPanel.setBorder(new EmptyBorder(10, 10, 10, 10));
        outerPanel.setLayout(new BorderLayout(0, 20));
        add(outerPanel);

        // Header message
        JPanel headerPanel = new JPanel();
        headerPanel.setLayout(new BoxLayout(headerPanel, BoxLayout.Y_AXIS));

        JLabel titleLabel = new JLabel("Track hubs");
        titleLabel.setFont(titleLabel.getFont().deriveFont(Font.BOLD, 16));
        headerPanel.add(titleLabel);
        headerPanel.add(Box.createVerticalStrut(20));

        JLabel headerMessage = new JLabel();
        headerMessage.setText("<html>Hubs selected from the UCSC Track Hub Registry or added by URL.  Check hubs to " +
                "remove and press Remove.<br>A removed hub is dropped from the 'Tracks' menu for every genome it " +
                "lists.  Hubs for the current genome (" + GenomeManager.getInstance().getCurrentGenome().getDisplayName() +
                ") are listed first.</html>");
        headerPanel.add(headerMessage);
        outerPanel.add(headerPanel, BorderLayout.NORTH);

        // Table
        HubTableModel model = new HubTableModel(this.hubDescriptors);
        HubTable table = new HubTable(model);

        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
        outerPanel.add(scrollPane, BorderLayout.CENTER);

        // Button panel
        JPanel buttonPanel = new JPanel();
        ((FlowLayout) buttonPanel.getLayout()).setAlignment(FlowLayout.RIGHT);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> {
            canceled = true;
            setVisible(false);
            dispose();
        });

        JButton removeButton = new JButton("Remove");
        removeButton.addActionListener(e -> {
            canceled = false;
            setVisible(false);
            dispose();
        });

        buttonPanel.add(removeButton);
        buttonPanel.add(cancelButton);
        outerPanel.add(buttonPanel, BorderLayout.SOUTH);
    }
}
