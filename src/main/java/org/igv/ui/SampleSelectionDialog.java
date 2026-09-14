package org.igv.ui;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Dialog for entering a list of sample IDs to show in a track.
 */
public class SampleSelectionDialog extends IGVDialog {

    private JTextArea textArea;
    private boolean isCanceled = true;

    /**
     * @param currentSamples IDs the text box starts with
     * @param allSamples     IDs restored by "Show All"
     */
    public SampleSelectionDialog(Frame parent, List<String> currentSamples, List<String> allSamples) {
        super(parent, "Filter Samples", true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        initComponents(currentSamples, allSamples);
        setLocationRelativeTo(parent);
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    /**
     * @return the entered sample IDs, or null if none were entered (show all samples)
     */
    public List<String> getSampleIds() {
        return parseSampleIds(textArea.getText());
    }

    /**
     * Parse sample IDs separated by newlines, tabs, or commas.  Spaces are not separators, as sample IDs may
     * contain them.
     */
    public static List<String> parseSampleIds(String text) {
        List<String> ids = Arrays.stream(text.split("[\\r\\n\\t,]+"))
                .map(String::trim)
                .filter(s -> !s.isEmpty())
                .distinct()
                .collect(Collectors.toList());
        return ids.isEmpty() ? null : ids;
    }

    private void initComponents(List<String> currentSamples, List<String> allSamples) {

        JPanel contentPanel = new JPanel(new BorderLayout(0, 8));
        contentPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        setContentPane(contentPanel);

        contentPanel.add(new JLabel("<html>Only the samples listed are shown, in list order.  Remove IDs to hide samples,<br>" +
                "or enter IDs one per line or separated by commas or tabs.  A filter by attribute also applies."),
                BorderLayout.NORTH);

        textArea = new JTextArea(15, 40);
        textArea.setText(String.join("\n", currentSamples));
        textArea.setCaretPosition(0);
        contentPanel.add(new JScrollPane(textArea), BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));

        JButton showAllButton = new JButton("Show All");
        showAllButton.addActionListener(evt -> {
            textArea.setText(String.join("\n", allSamples));
            textArea.setCaretPosition(0);
        });
        buttonPanel.add(showAllButton);

        JButton okButton = new JButton("OK");
        okButton.addActionListener(evt -> {
            isCanceled = false;
            setVisible(false);
            dispose();
        });
        buttonPanel.add(okButton);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(evt -> {
            isCanceled = true;
            setVisible(false);
            dispose();
        });
        buttonPanel.add(cancelButton);

        contentPanel.add(buttonPanel, BorderLayout.SOUTH);
        getRootPane().setDefaultButton(okButton);

        pack();
    }
}
