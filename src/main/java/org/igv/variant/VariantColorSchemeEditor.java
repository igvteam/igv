package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ContinuousColorScale;
import org.igv.ui.IGV;
import org.igv.ui.color.ColorSwatch;
import org.igv.ui.legend.HeatmapLegendEditor;
import org.igv.ui.util.MessageUtils;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Window;
import java.util.Map;

/**
 * Edits a saved color scheme -- every attribute it covers, not just the one a track happens to be colored by.
 * Reached from the Variants tab of the preferences editor, so a scheme can be adjusted without loading a file
 * that uses it.
 * <p>
 * Editing a scheme shipped with IGV saves a copy the user owns, which then shadows the original.
 */
public class VariantColorSchemeEditor extends JDialog {

    private static Logger log = LogManager.getLogger(VariantColorSchemeEditor.class);

    private final VariantColorScheme scheme;
    private final JPanel contentPanel = new JPanel();
    private boolean saved = false;

    public VariantColorSchemeEditor(Window owner, VariantColorScheme scheme) {

        super(owner, "Edit " + scheme.getName(), ModalityType.APPLICATION_MODAL);
        this.scheme = scheme;

        JPanel content = new JPanel(new BorderLayout(0, 8));
        content.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        JLabel description = new JLabel(scheme.isBuiltIn()
                ? "<html>This scheme ships with IGV.  Saving writes your own copy, which takes precedence over it."
                : "<html>Colors for the values of VCF INFO attributes.");
        content.add(description, BorderLayout.NORTH);

        contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
        JScrollPane scrollPane = new JScrollPane(contentPanel);
        scrollPane.setPreferredSize(new Dimension(420, 320));
        scrollPane.getVerticalScrollBar().setUnitIncrement(16);
        content.add(scrollPane, BorderLayout.CENTER);

        JButton saveButton = new JButton("Save");
        saveButton.addActionListener(e -> save());

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> dispose());

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 0));
        buttonPanel.add(saveButton);
        buttonPanel.add(cancelButton);
        content.add(buttonPanel, BorderLayout.SOUTH);

        setContentPane(content);
        populate();
        pack();
        setLocationRelativeTo(owner);
    }

    /**
     * @return true if the scheme was saved, so the caller can refresh.
     */
    public boolean isSaved() {
        return saved;
    }

    private void populate() {

        contentPanel.removeAll();

        for (String key : scheme.getKeys()) {

            JLabel heading = new JLabel(key);
            heading.setFont(heading.getFont().deriveFont(Font.BOLD));
            heading.setBorder(BorderFactory.createEmptyBorder(8, 2, 2, 2));
            heading.setAlignmentX(LEFT_ALIGNMENT);
            contentPanel.add(heading);

            AbstractColorScale scale = scheme.getScale(key);
            if (scale instanceof ContinuousColorScale) {
                contentPanel.add(createScalePanel(key, (ContinuousColorScale) scale));
                continue;
            }

            if (scheme.isCategorical(key)) {
                contentPanel.add(note("Values are categories, not a range."));
            }

            for (Map.Entry<String, Color> entry : scheme.getColors(key).entrySet()) {
                contentPanel.add(createColorRow(key, entry.getKey(), entry.getValue()));
            }

            Color defaultColor = scheme.getDefaultColor(key);
            if (defaultColor != null) {
                contentPanel.add(createColorRow(key, VariantColorScheme.WILDCARD, defaultColor));
            }

            if (scheme.getColors(key).isEmpty() && defaultColor == null) {
                contentPanel.add(note("No colors set -- values take them from the palette."));
            }
        }

        contentPanel.add(Box.createVerticalGlue());
        contentPanel.revalidate();
        contentPanel.repaint();
    }

    private JPanel createColorRow(String key, String value, Color color) {

        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));
        row.setAlignmentX(LEFT_ALIGNMENT);

        ColorSwatch swatch = new ColorSwatch(color);
        swatch.addColorChangeListener(c -> scheme.setColor(key, value, c));
        row.add(swatch);
        row.add(new JLabel(VariantColorScheme.WILDCARD.equals(value) ? "* (all other values)" : value));

        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, row.getPreferredSize().height));
        return row;
    }

    private JPanel createScalePanel(String key, ContinuousColorScale scale) {

        JPanel panel = new JPanel(new BorderLayout(6, 4));
        panel.setBorder(BorderFactory.createEmptyBorder(2, 4, 4, 4));
        panel.setAlignmentX(LEFT_ALIGNMENT);

        JPanel gradient = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                ContinuousColorScale current = (ContinuousColorScale) scheme.getScale(key);
                double min = current.getMinimum();
                double max = current.getMaximum();
                for (int x = 0; x < getWidth(); x++) {
                    double value = min + (max - min) * x / Math.max(1, getWidth() - 1);
                    g.setColor(current.getColor((float) value));
                    g.drawLine(x, 0, x, getHeight());
                }
            }
        };
        gradient.setPreferredSize(new Dimension(260, 22));
        panel.add(gradient, BorderLayout.CENTER);

        JButton edit = new JButton("Edit Scale...");
        edit.addActionListener(e -> {
            ContinuousColorScale current = (ContinuousColorScale) scheme.getScale(key);
            HeatmapLegendEditor editor = new HeatmapLegendEditor(
                    getOwner() instanceof java.awt.Frame ? (java.awt.Frame) getOwner() : null, true, current);
            editor.setTitle("Color scale for " + key);
            editor.setVisible(true);
            if (!editor.isCanceled()) {
                scheme.setScale(key, editor.getColorScheme());
                gradient.repaint();
            }
        });
        panel.add(edit, BorderLayout.EAST);

        panel.setMaximumSize(new Dimension(Integer.MAX_VALUE, panel.getPreferredSize().height));
        return panel;
    }

    private JLabel note(String text) {
        JLabel label = new JLabel(text);
        label.setFont(label.getFont().deriveFont(Font.ITALIC, label.getFont().getSize() - 1f));
        label.setBorder(BorderFactory.createEmptyBorder(2, 6, 2, 2));
        label.setAlignmentX(LEFT_ALIGNMENT);
        return label;
    }

    private void save() {
        try {
            if (scheme.isBuiltIn()) {
                // Cannot write to a scheme on the classpath -- save the user their own copy
                scheme.setName(scheme.getName() + " (edited)");
            }
            VariantColorSchemes.save(scheme);
            saved = true;
            if (IGV.hasInstance()) {
                IGV.getInstance().repaint();
            }
            dispose();
        } catch (Exception e) {
            log.error("Error saving color scheme " + scheme.getName(), e);
            MessageUtils.showMessage("Error saving color scheme: " + e.getMessage());
        }
    }
}
