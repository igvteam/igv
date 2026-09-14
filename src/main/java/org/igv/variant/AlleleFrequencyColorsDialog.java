package org.igv.variant;

import org.igv.renderer.ColorStopScale;
import org.igv.ui.IGV;
import org.igv.ui.color.ColorSwatch;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

/**
 * Edits the color stops used to color variants by allele frequency, and serves as their legend.  Frequencies are
 * shown as percentages.  Changes apply to all variant tracks and are saved as a preference.
 */
public class AlleleFrequencyColorsDialog extends JDialog {

    private final List<StopRow> rows = new ArrayList<>();
    private final JPanel stopsPanel = new JPanel();
    private final ScalePreview preview = new ScalePreview();

    private static class StopRow {
        JTextField percent;
        Color color;
    }

    public AlleleFrequencyColorsDialog(Frame owner) {
        super(owner, "Allele Frequency Colors", true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        initComponents();
        load(AlleleFrequencyColors.getScale());
        pack();
        setLocationRelativeTo(owner);
    }

    private void initComponents() {

        JPanel content = new JPanel(new BorderLayout(0, 8));
        content.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        setContentPane(content);

        content.add(new JLabel("<html>Variants are colored by allele frequency, blended on a log scale between stops.<br>" +
                "Variants with no frequency take the color of the lowest stop."), BorderLayout.NORTH);

        stopsPanel.setLayout(new BoxLayout(stopsPanel, BoxLayout.Y_AXIS));
        JPanel center = new JPanel(new BorderLayout(0, 8));
        center.add(preview, BorderLayout.NORTH);
        center.add(stopsPanel, BorderLayout.CENTER);
        content.add(center, BorderLayout.CENTER);

        JButton addButton = new JButton("Add Stop");
        addButton.addActionListener(e -> addStop());
        JButton resetButton = new JButton("Reset to Defaults");
        resetButton.addActionListener(e -> load(AlleleFrequencyColors.DEFAULT_SCALE));
        JPanel leftButtons = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        leftButtons.add(addButton);
        leftButtons.add(resetButton);

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> save());
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> dispose());
        JPanel rightButtons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 4, 0));
        rightButtons.add(okButton);
        rightButtons.add(cancelButton);

        JPanel buttons = new JPanel(new BorderLayout());
        buttons.add(leftButtons, BorderLayout.WEST);
        buttons.add(rightButtons, BorderLayout.EAST);
        content.add(buttons, BorderLayout.SOUTH);

        getRootPane().setDefaultButton(okButton);
    }

    /**
     * Show a scale's stops, highest frequency first.
     */
    private void load(ColorStopScale scale) {
        rows.clear();
        List<ColorStopScale.Stop> stops = scale.getStops();
        for (int i = stops.size() - 1; i >= 0; i--) {
            addRow(stops.get(i).value(), stops.get(i).color());
        }
        refresh();
    }

    /**
     * Add a stop a factor of 10 below the lowest, in its color.
     */
    private void addStop() {
        ColorStopScale scale = currentScale();
        ColorStopScale.Stop lowest = scale == null ? null : scale.getStops().get(0);
        addRow(lowest == null ? 0.01 : lowest.value() / 10, lowest == null ? Color.GRAY : lowest.color());
        refresh();
    }

    private void addRow(double frequency, Color color) {
        StopRow row = new StopRow();
        row.percent = new JTextField(toPercent(frequency), 8);
        row.percent.setHorizontalAlignment(JTextField.RIGHT);
        row.percent.getDocument().addDocumentListener(new DocumentListener() {
            public void insertUpdate(DocumentEvent e) { updatePreview(); }
            public void removeUpdate(DocumentEvent e) { updatePreview(); }
            public void changedUpdate(DocumentEvent e) { updatePreview(); }
        });
        row.color = color;
        rows.add(row);
    }

    private void refresh() {

        stopsPanel.removeAll();
        for (StopRow row : rows) {

            ColorSwatch swatch = new ColorSwatch(row.color);
            swatch.addColorChangeListener(c -> {
                row.color = c;
                updatePreview();
            });

            JButton remove = new JButton("×");
            remove.setToolTipText("Remove this stop");
            remove.setMargin(new Insets(0, 0, 0, 0));
            remove.setPreferredSize(new Dimension(22, 20));
            remove.setEnabled(rows.size() > 1);
            remove.addActionListener(e -> {
                rows.remove(row);
                refresh();
            });

            JPanel left = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));
            left.add(swatch);
            left.add(row.percent);
            left.add(new JLabel("%"));

            JPanel right = new JPanel(new FlowLayout(FlowLayout.RIGHT, 4, 2));
            right.add(remove);

            JPanel line = new JPanel(new BorderLayout());
            line.add(left, BorderLayout.WEST);
            line.add(right, BorderLayout.EAST);
            line.setMaximumSize(new Dimension(Integer.MAX_VALUE, line.getPreferredSize().height));
            stopsPanel.add(line);
        }
        stopsPanel.revalidate();
        updatePreview();
        if (isDisplayable()) {
            pack();
        }
    }

    private void updatePreview() {
        preview.setScale(currentScale());
    }

    /**
     * @return the scale the rows describe, or null if a percentage isn't a positive number
     */
    private ColorStopScale currentScale() {
        List<ColorStopScale.Stop> stops = new ArrayList<>();
        for (StopRow row : rows) {
            Double frequency = fromPercent(row.percent.getText());
            if (frequency == null) {
                return null;
            }
            stops.add(new ColorStopScale.Stop(frequency, row.color));
        }
        return stops.isEmpty() ? null : new ColorStopScale(stops);
    }

    private void save() {
        ColorStopScale scale = currentScale();
        if (scale == null) {
            JOptionPane.showMessageDialog(this, "Each stop needs a frequency greater than 0%.",
                    "Allele Frequency Colors", JOptionPane.WARNING_MESSAGE);
            return;
        }
        AlleleFrequencyColors.setScale(scale);
        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
        dispose();
    }

    static String toPercent(double frequency) {
        return BigDecimal.valueOf(frequency).movePointRight(2).stripTrailingZeros().toPlainString();
    }

    /**
     * @return the frequency for a percentage, or null if it isn't a positive number
     */
    static Double fromPercent(String text) {
        try {
            double frequency = new BigDecimal(text.trim()).movePointLeft(2).doubleValue();
            return frequency > 0 ? frequency : null;
        } catch (NumberFormatException e) {
            return null;
        }
    }

    /**
     * The scale as a bar, from the highest stop down to a factor of 10 below the lowest, with each stop marked.
     */
    private static class ScalePreview extends JComponent {

        private static final int BAR_HEIGHT = 20;
        private ColorStopScale scale;

        ScalePreview() {
            setPreferredSize(new Dimension(400, BAR_HEIGHT + 22));
        }

        void setScale(ColorStopScale scale) {
            this.scale = scale;
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g) {

            Graphics2D g2 = (Graphics2D) g.create();
            int x0 = 10;
            int width = getWidth() - 20;
            if (scale == null || width <= 0) {
                g2.setColor(Color.GRAY);
                g2.drawRect(x0, 0, Math.max(width, 1), BAR_HEIGHT);
                g2.dispose();
                return;
            }

            List<ColorStopScale.Stop> stops = scale.getStops();
            double logHigh = Math.log10(Math.max(stops.get(stops.size() - 1).value(), stops.get(0).value() * 10));
            double logLow = Math.log10(stops.get(0).value() / 10);

            for (int i = 0; i < width; i++) {
                double frequency = Math.pow(10, logHigh - (logHigh - logLow) * i / width);
                g2.setColor(scale.getColor(frequency));
                g2.fillRect(x0 + i, 0, 1, BAR_HEIGHT);
            }

            g2.setColor(getForeground());
            g2.setFont(g2.getFont().deriveFont(g2.getFont().getSize() - 2f));
            FontMetrics fm = g2.getFontMetrics();
            int lastLabelEnd = Integer.MIN_VALUE;
            for (int i = stops.size() - 1; i >= 0; i--) {
                double position = (logHigh - Math.log10(stops.get(i).value())) / (logHigh - logLow);
                int x = x0 + (int) Math.round(position * width);
                g2.drawLine(x, BAR_HEIGHT, x, BAR_HEIGHT + 3);
                String label = toPercent(stops.get(i).value()) + "%";
                int labelX = Math.max(0, Math.min(getWidth() - fm.stringWidth(label), x - fm.stringWidth(label) / 2));
                if (labelX > lastLabelEnd + 4) {        // Skip labels that would overlap
                    g2.drawString(label, labelX, BAR_HEIGHT + 4 + fm.getAscent());
                    lastLabelEnd = labelX + fm.stringWidth(label);
                }
            }
            g2.dispose();
        }
    }
}
