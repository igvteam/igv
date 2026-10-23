package org.igv.ucsc.hub;

import org.igv.feature.genome.load.GenomeConfig;
import org.igv.ui.genome.GenomeListItem;
import org.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Offers the genomes referenced by a multi-genome track hub, that are not already in the genome list, for addition
 * to the list.  All are checked by default.  Ids the hub references that IGV does not know are listed in a note.
 */
public class HubGenomeSelectionDialog extends JDialog {

    /**
     * A genome referenced by the hub.  Exactly one of config (an assembly defined by the hub) or item (a hosted genome
     * record) is set.
     */
    public record Candidate(String id, String name, GenomeConfig config, GenomeListItem item) {
    }

    private final List<Candidate> candidates;
    private final List<JCheckBox> checkBoxes = new ArrayList<>();
    private boolean canceled = true;

    public HubGenomeSelectionDialog(Frame owner, Hub hub, List<Candidate> candidates, List<String> unknownIds) {

        super(owner, "Add Genomes from Hub", true);
        this.candidates = candidates;

        JPanel outerPanel = new JPanel(new BorderLayout(0, 10));
        outerPanel.setBorder(new EmptyBorder(10, 10, 10, 10));
        add(outerPanel);

        // Header: hub name and prompt
        JPanel headerPanel = new JPanel();
        headerPanel.setLayout(new BoxLayout(headerPanel, BoxLayout.Y_AXIS));
        headerPanel.setBorder(new EmptyBorder(0, 0, 5, 0));
        JLabel titleLabel = new JLabel(hub.getShortLabel());
        titleLabel.setFont(titleLabel.getFont().deriveFont(Font.BOLD, 14));
        titleLabel.setAlignmentX(LEFT_ALIGNMENT);
        headerPanel.add(titleLabel);
        if (!candidates.isEmpty()) {
            headerPanel.add(Box.createVerticalStrut(8));
            JLabel message = new JLabel("Add the genomes referenced by this hub to the genome selection list?");
            message.setAlignmentX(LEFT_ALIGNMENT);
            headerPanel.add(message);
        }
        outerPanel.add(headerPanel, BorderLayout.NORTH);

        // Selection list: select all / none above the scrolling list of genomes
        JPanel selectionPanel = new JPanel(new BorderLayout(0, 5));
        selectionPanel.setBorder(BorderFactory.createTitledBorder("Genomes"));

        JPanel listPanel = new JPanel();
        listPanel.setLayout(new BoxLayout(listPanel, BoxLayout.Y_AXIS));
        listPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        listPanel.setBackground(UIManager.getColor("List.background"));
        for (Candidate candidate : candidates) {
            String label = candidate.name().contains(candidate.id()) ?
                    candidate.name() : candidate.name() + " (" + candidate.id() + ")";
            JCheckBox cb = new JCheckBox(label, true);
            cb.setOpaque(false);
            checkBoxes.add(cb);
            listPanel.add(cb);
        }

        if (!candidates.isEmpty()) {
            JPanel selectPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
            JCheckBox selectAll = new JCheckBox("Select all", true);
            selectAll.addActionListener(e -> checkBoxes.forEach(cb -> cb.setSelected(selectAll.isSelected())));
            selectPanel.add(selectAll);
            selectionPanel.add(selectPanel, BorderLayout.NORTH);
        }

        JScrollPane scrollPane = new JScrollPane(listPanel);
        scrollPane.setBorder(BorderFactory.createLineBorder(UIManager.getColor("Separator.foreground")));
        scrollPane.getViewport().setBackground(UIManager.getColor("List.background"));
        selectionPanel.add(scrollPane, BorderLayout.CENTER);
        outerPanel.add(selectionPanel, BorderLayout.CENTER);

        // Footer: note on ids IGV does not know, and the buttons
        JPanel footerPanel = new JPanel(new BorderLayout(0, 5));
        footerPanel.setBorder(new EmptyBorder(5, 0, 0, 0));

        if (!unknownIds.isEmpty()) {
            JLabel note = new JLabel("<html><i>Not known to IGV and cannot be added: " + String.join(", ", unknownIds) + "</i></html>");
            footerPanel.add(note, BorderLayout.CENTER);
        }

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 0));
        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> {
            canceled = false;
            setVisible(false);
        });
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> {
            canceled = true;
            setVisible(false);
        });
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        getRootPane().setDefaultButton(okButton);
        footerPanel.add(buttonPanel, BorderLayout.SOUTH);
        outerPanel.add(footerPanel, BorderLayout.SOUTH);

        setSize(500, Math.min(600, 230 + 25 * candidates.size()));
        if (owner != null) {
            setLocationRelativeTo(owner);
        }
    }

    public boolean isCanceled() {
        return canceled;
    }

    public List<Candidate> getSelectedCandidates() {
        List<Candidate> selected = new ArrayList<>();
        for (int i = 0; i < candidates.size(); i++) {
            if (checkBoxes.get(i).isSelected()) {
                selected.add(candidates.get(i));
            }
        }
        return selected;
    }

    /**
     * Offer a switch to one of the genomes referenced by a hub, with a warning that switching clears the loaded
     * tracks.  Returns the chosen genome, or null if the user cancels.
     */
    public static GenomeListItem selectGenomeToSwitchTo(Frame owner, List<GenomeListItem> hubGenomes) {
        String message = "The current genome is not referenced by this hub.  Switch to one of the hub's genomes?\n" +
                "Switching genomes will clear the currently loaded tracks.";
        Object[] options = hubGenomes.toArray();
        final Object[] choice = new Object[1];
        UIUtilities.invokeAndWaitOnEventThread(() -> choice[0] = JOptionPane.showInputDialog(owner, message,
                "Switch Genome", JOptionPane.QUESTION_MESSAGE, null, options, options[0]));
        return (GenomeListItem) choice[0];
    }
}
