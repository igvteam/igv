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
import java.awt.Frame;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Shows the color assigned to each value of the INFO attribute a variant track is colored by, and lets the user
 * change them.
 * <p>
 * This is where colors are actually chosen, rather than in preferences: the attribute's values are not knowable
 * in advance, but they are knowable here -- they are harvested from the loaded variants, which is also the point
 * at which the user can see what needs recoloring.  Changes apply to this track only; "Save as Scheme" promotes
 * them to a scheme that applies to every VCF with this attribute.
 * <p>
 * The dialog is modeless so the track repaints under it as colors change.
 */
public class VariantColorLegendDialog extends JDialog {

    private static Logger log = LogManager.getLogger(VariantColorLegendDialog.class);

    private final VariantTrack track;
    private final String infoKey;
    private final JPanel valuePanel = new JPanel();
    private final JLabel description = new JLabel();
    private final JButton editScaleButton = new JButton("Edit Scale...");
    private final JButton saveButton = new JButton("Save as Scheme...");
    private final JButton resetButton = new JButton("Reset");

    public VariantColorLegendDialog(Frame owner, VariantTrack track, String infoKey) {

        super(owner, "Colors for " + infoKey, false);
        this.track = track;
        this.infoKey = infoKey;

        JPanel content = new JPanel(new BorderLayout(0, 8));
        content.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        content.add(description, BorderLayout.NORTH);

        valuePanel.setLayout(new BoxLayout(valuePanel, BoxLayout.Y_AXIS));
        JScrollPane scrollPane = new JScrollPane(valuePanel);
        scrollPane.setPreferredSize(new Dimension(380, 260));
        scrollPane.getVerticalScrollBar().setUnitIncrement(16);
        content.add(scrollPane, BorderLayout.CENTER);

        editScaleButton.setToolTipText("Change the range and colors of the scale");
        editScaleButton.addActionListener(e -> editScale());

        saveButton.setToolTipText("Save these colors so they apply to every VCF with this attribute");
        saveButton.addActionListener(e -> saveAsScheme());

        resetButton.setToolTipText("Discard the colors chosen for this track");
        resetButton.addActionListener(e -> {
            track.clearAttributeColorOverrides(infoKey);
            repaintTrack();
            populate();
        });

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> dispose());

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 0));
        buttonPanel.add(editScaleButton);
        buttonPanel.add(saveButton);
        buttonPanel.add(resetButton);
        buttonPanel.add(closeButton);
        content.add(buttonPanel, BorderLayout.SOUTH);

        setContentPane(content);
        populate();
        pack();
        setLocationRelativeTo(owner);
    }

    /**
     * Build a row per value -- see {@link #valuesToShow()}.
     */
    private void populate() {

        valuePanel.removeAll();

        AbstractColorScale scale = VariantColorSchemes.getScale(infoKey);

        description.setText("<html>Values of <b>" + infoKey + "</b> in the current view.  "
                + (scale == null ? "Click a color to change it for this track." : "Colored by a scale."));

        editScaleButton.setVisible(scale instanceof ContinuousColorScale);
        saveButton.setVisible(scale == null);
        resetButton.setVisible(scale == null);

        if (scale != null) {
            valuePanel.add(createScaleRow(scale));
            valuePanel.add(Box.createVerticalGlue());
            valuePanel.revalidate();
            valuePanel.repaint();
            return;
        }

        Set<String> values = valuesToShow();

        if (values.isEmpty()) {
            JLabel empty = new JLabel("No variants in view have a value for " + infoKey + ".");
            empty.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
            valuePanel.add(empty);
        } else {
            for (String value : values) {
                valuePanel.add(createRow(value));
            }
        }

        valuePanel.add(Box.createVerticalGlue());
        valuePanel.revalidate();
        valuePanel.repaint();
    }

    /**
     * The values to list, and to write on "Save as Scheme": those in view, plus any the user chose a color for
     * on this track, which may have scrolled out of view.  Sorted case insensitively, keeping the spelling seen
     * in view.
     * <p>
     * Colors IGV assigned from the palette to values that have scrolled out of view are deliberately left out.
     * That table grows with every value ever drawn -- thousands, for an attribute like a gene name -- and keys
     * it by lower case value.  The assignments are not lost: the track keeps them, they are saved in sessions,
     * and they come back when the values do.
     */
    private Set<String> valuesToShow() {
        Set<String> values = new TreeSet<>(String.CASE_INSENSITIVE_ORDER);
        values.addAll(track.getAttributeValues(infoKey));
        values.addAll(track.getAttributeColorOverrides(infoKey).keySet());
        return values;
    }

    /**
     * A numeric attribute is colored by a scale rather than by value, so the legend is the gradient itself.
     */
    private JPanel createScaleRow(AbstractColorScale scale) {

        JPanel panel = new JPanel(new BorderLayout(0, 4));
        panel.setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));

        double[] range = track.getAttributeRange(infoKey);
        JLabel label = new JLabel(range == null ? "Colored by a scale" :
                String.format("Values in view: %s to %s",
                        ColorScaleBar.format(range[0], range[1] - range[0]),
                        ColorScaleBar.format(range[1], range[1] - range[0])));
        panel.add(label, BorderLayout.NORTH);

        panel.add(new ColorScaleBar(scale instanceof ContinuousColorScale ? (ContinuousColorScale) scale : null, 340),
                BorderLayout.CENTER);

        panel.setMaximumSize(new Dimension(Integer.MAX_VALUE, panel.getPreferredSize().height));
        return panel;
    }

    private JPanel createRow(String value) {

        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));

        ColorSwatch swatch = new ColorSwatch(track.getAttributeColor(infoKey, value));
        swatch.addColorChangeListener(color -> {
            track.setAttributeColorOverride(infoKey, value, color);
            repaintTrack();
        });
        row.add(swatch);

        row.add(new JLabel(value));

        // Say where the color comes from, so it is clear what "Reset" and "Save as Scheme" will do
        String origin = colorOrigin(value);
        if (origin != null) {
            JLabel originLabel = new JLabel(origin);
            originLabel.setFont(originLabel.getFont().deriveFont(Font.ITALIC, originLabel.getFont().getSize() - 1f));
            row.add(originLabel);
        }

        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, row.getPreferredSize().height));
        return row;
    }

    private String colorOrigin(String value) {
        if (track.getAttributeColorOverrides(infoKey).containsKey(value)) {
            return "(this track)";
        }
        for (VariantColorScheme scheme : VariantColorSchemes.getSchemes()) {
            if (scheme.getColor(infoKey, value) != null) {
                return "(" + scheme.getName() + ")";
            }
        }
        return "(assigned)";
    }

    /**
     * Reopen the scale editor and save the result back over the scheme it came from, so the change applies
     * wherever the scale did.
     */
    private void editScale() {

        AbstractColorScale scale = VariantColorSchemes.getScale(infoKey);
        if (!(scale instanceof ContinuousColorScale)) {
            return;
        }

        ContinuousColorScale edited =
                HeatmapLegendEditor.edit(this, "Color scale for " + infoKey, (ContinuousColorScale) scale);
        if (edited == null) {
            return;
        }

        VariantColorScheme existing = VariantColorSchemes.getSchemeForScale(infoKey);
        VariantColorScheme updated = existing == null ? new VariantColorScheme(infoKey + " colorscale") : existing.copy();
        if (updated.isBuiltIn()) {
            updated.setName(infoKey + " colorscale");
        }
        updated.setScale(infoKey, edited);

        try {
            VariantColorSchemes.save(updated);
            repaintTrack();
            populate();
        } catch (Exception e) {
            log.error("Error saving color scale for " + infoKey, e);
            MessageUtils.showMessage("Error saving color scale: " + e.getMessage());
        }
    }

    private void saveAsScheme() {

        String name = MessageUtils.showInputDialog("Scheme name", infoKey + " colors");
        if (name == null || name.trim().isEmpty()) {
            return;
        }

        Map<String, Color> colors = new LinkedHashMap<>();
        for (String value : valuesToShow()) {
            colors.put(value, track.getAttributeColor(infoKey, value));
        }

        try {
            VariantColorSchemes.saveScheme(name.trim(), infoKey, colors);
            // The colors now come from the scheme, so drop the per-track copies -- otherwise editing the scheme
            // later would have no visible effect on this track.
            track.clearAttributeColorOverrides(infoKey);
            repaintTrack();
            populate();
        } catch (Exception e) {
            log.error("Error saving variant color scheme", e);
            MessageUtils.showMessage("Error saving color scheme: " + e.getMessage());
        }
    }

    private void repaintTrack() {
        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
    }
}
