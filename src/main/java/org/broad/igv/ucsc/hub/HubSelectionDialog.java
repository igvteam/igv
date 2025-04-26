package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.ui.commandbar.GenomeListManager;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class HubSelectionDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    private static Map<String, HubSelectionDialog> instances = new java.util.HashMap<>();

    List<HubDescriptor> hubDescriptors;
    private boolean canceled;
    private JTable table;

    public HubSelectionDialog(Frame owner) {

        super(owner, "Select Track Hubs", true);

        setSize(owner.getWidth() - 20, 600);
        if (owner != null) {
            setLocationRelativeTo(owner);
        }


        List<HubDescriptor> allHubDescriptors = HubRegistry.getAllHubDescriptors();
        final List<HubDescriptor> allSelectedHubs = HubRegistry.getAllSelectedHubs();

        // Filter track hubs to those supporting current genome list or selected
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
    }

    private List<HubDescriptor> filterHubs(List<HubDescriptor> hubDescriptors, List<HubDescriptor> allSelectedHubs) {

        Set<String> allSelectedURLs = allSelectedHubs.stream().map(HubDescriptor::getUrl).collect(Collectors.toSet());

        return hubDescriptors.stream().filter(h -> {
            Set<String> genomeIDs = getGenomeIDList();
            String[] hubGenomes = h.getDbList().split(",");
            if (hubGenomes == null || hubGenomes.length == 0) {
                return true;
            }
            if(allSelectedURLs.contains(h.getUrl())) {
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
    private Set<String> getGenomeIDList() {
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
        return IntStream.range(0, table.getRowCount())
                .filter(i -> (Boolean) table.getValueAt(i, 0))
                .mapToObj(hubDescriptors::get)
                .collect(Collectors.toList());
    }


    private void initComponents(Collection<HubDescriptor> loadedHubs) {

        this.setLayout(new BorderLayout());

        Set<String> loadedURLs = loadedHubs.stream().map(HubDescriptor::getUrl).collect(Collectors.toSet());

        String[] columnNames = {"", "Name", "Description", "DB List", "Info"};
        Object[][] data = new Object[this.hubDescriptors.size()][5];
        for (int i = 0; i < this.hubDescriptors.size(); i++) {
            HubDescriptor hub = this.hubDescriptors.get(i);
            data[i][0] = loadedURLs.contains(hub.getUrl()); // Default unchecked checkbox
            data[i][1] = hub.getShortLabel();
            data[i][2] = hub.getLongLabel();
            data[i][3] = hub.getDbList();
            data[i][4] = hub.getDescriptionUrl();
        }

        table = new HubTable(data, columnNames);

        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setBorder(BorderFactory.createEmptyBorder(5,10,5,10));
        add(scrollPane, BorderLayout.CENTER);

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
        add(buttonPanel, BorderLayout.SOUTH);
    }
}

