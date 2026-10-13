package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ColorStopScale;
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
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Shows the color assigned to each value of the INFO attribute a variant track is colored by, and lets the user
 * change them.
 * <p>
 * This is where colors are actually chosen, rather than in preferences: the attribute's values are not knowable
 * in advance, but they are knowable here -- they are harvested from the loaded variants, which is also the point
 * at which the user can see what needs recoloring.  Nothing changes until "OK", which saves the colors as a
 * scheme that applies to every VCF with this attribute; "Cancel" discards them.
 */
public class VariantColorLegendDialog extends JDialog {

    private static Logger log = LogManager.getLogger(VariantColorLegendDialog.class);

    private final VariantTrack track;
    private final String infoKey;
    private final JPanel valuePanel = new JPanel();
    private final JLabel description = new JLabel();
    private final JButton editScaleButton = new JButton("Edit Scale...");

    // Edits not yet saved -- applied on "OK", discarded on "Cancel"
    private final Map<String, Color> pendingColors = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);
    private ContinuousColorScale pendingScale;

    public VariantColorLegendDialog(Frame owner, VariantTrack track, String infoKey) {

        super(owner, "Colors for " + infoKey, true);
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

        JButton okButton = new JButton("OK");
        okButton.setToolTipText("Save these colors so they apply to every VCF with this attribute");
        okButton.addActionListener(e -> {
            if (save()) {
                dispose();
            }
        });

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> dispose());

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 0));
        buttonPanel.add(editScaleButton);
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
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

        AbstractColorScale scale = pendingScale != null ? pendingScale : VariantColorSchemes.getScale(infoKey);

        description.setText("<html>Values of <b>" + infoKey + "</b> among the loaded features.  "
                + (scale == null ? "Click a color to change it." : "Colored by a scale."));

        editScaleButton.setVisible(scale instanceof ContinuousColorScale || scale instanceof ColorStopScale);

        if (scale != null) {
            valuePanel.add(createScaleRow(scale));
            valuePanel.add(Box.createVerticalGlue());
            valuePanel.revalidate();
            valuePanel.repaint();
            return;
        }

        Set<String> values = valuesToShow();

        if (values.isEmpty()) {
            JLabel empty = new JLabel("No loaded features have a value for " + infoKey + ".");
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
     * The values to list, and to write on "OK": those among the loaded features (see
     * {@link VariantTrack#getAttributeValues}), plus any the user chose a color for on this track, which may no
     * longer be loaded.  Sorted case insensitively, keeping the spelling seen in the loaded features.
     * <p>
     * Colors IGV assigned from the palette to values that are no longer loaded are deliberately left out.
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
                String.format("Values in loaded features: %s to %s",
                        ColorScaleBar.format(range[0], range[1] - range[0]),
                        ColorScaleBar.format(range[1], range[1] - range[0])));
        panel.add(label, BorderLayout.NORTH);

        if (scale instanceof ColorStopScale) {
            // A well-known allele frequency field, colored by rarity
            AlleleFrequencyColorsDialog.ScalePreview preview = new AlleleFrequencyColorsDialog.ScalePreview();
            preview.setScale((ColorStopScale) scale);
            panel.add(preview, BorderLayout.CENTER);
        } else {
            panel.add(new ColorScaleBar(scale instanceof ContinuousColorScale ? (ContinuousColorScale) scale : null, 340),
                    BorderLayout.CENTER);
        }

        panel.setMaximumSize(new Dimension(Integer.MAX_VALUE, panel.getPreferredSize().height));
        return panel;
    }

    private JPanel createRow(String value) {

        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));

        ColorSwatch swatch = new ColorSwatch(currentColor(value));
        swatch.addColorChangeListener(color -> pendingColors.put(value, color));
        row.add(swatch);

        row.add(new JLabel(value));

        // Say where the color comes from, so it is clear what "OK" will change
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

    private Color currentColor(String value) {
        Color pending = pendingColors.get(value);
        return pending != null ? pending : track.getAttributeColor(infoKey, value);
    }

    /**
     * Reopen the scale editor.  The result is kept until "OK".
     */
    private void editScale() {

        AbstractColorScale scale = pendingScale != null ? pendingScale : VariantColorSchemes.getScale(infoKey);
        if (scale instanceof ColorStopScale) {
            // Rarity colors are shared by all allele frequency coloring, and edited in one place
            new AlleleFrequencyColorsDialog(getOwner() instanceof Frame frame ? frame : null).setVisible(true);
            populate();
            repaintTrack();
            return;
        }
        if (!(scale instanceof ContinuousColorScale)) {
            return;
        }

        ContinuousColorScale edited =
                HeatmapLegendEditor.edit(this, "Color scale for " + infoKey, (ContinuousColorScale) scale);
        if (edited != null) {
            pendingScale = edited;
            populate();
        }
    }

    /**
     * Save the edits back over the scheme the colors came from, so the change applies wherever that scheme did.
     * A scheme shipped with IGV cannot be written, so it is copied under a default name, which then shadows it.
     * Updating the existing scheme rather than writing a new one keeps what it holds beyond these edits -- colors
     * for values not loaded now, a categorical declaration.
     *
     * @return true if saved, or there was nothing to save
     */
    private boolean save() {

        if (pendingScale == null && pendingColors.isEmpty()) {
            return true;
        }

        VariantColorScheme existing = pendingScale != null
                ? VariantColorSchemes.getSchemeForScale(infoKey) : getSchemeForColors();
        String defaultName = infoKey + (pendingScale != null ? " colorscale" : " colors");
        VariantColorScheme updated = existing == null ? new VariantColorScheme(defaultName) : existing.copy();
        if (updated.isBuiltIn()) {
            updated.setName(defaultName);
        }

        if (pendingScale != null) {
            updated.setScale(infoKey, pendingScale);
        } else {
            // All listed values, not just the edited ones -- colors assigned from the palette depend on the order
            // values are encountered, so they are fixed here along with the user's choices
            for (String value : valuesToShow()) {
                updated.setColor(infoKey, value, currentColor(value));
            }
            if (track.isNumericAttribute(infoKey)) {
                // Colored by value, so it holds categories -- without the declaration it would be asked for a scale
                updated.setCategorical(infoKey);
            }
        }

        try {
            VariantColorSchemes.save(updated);
            // The colors now come from the scheme, so drop any per-track copies (from a session) -- otherwise they
            // would hide the change.
            track.clearAttributeColorOverrides(infoKey);
            repaintTrack();
            return true;
        } catch (Exception e) {
            log.error("Error saving colors for " + infoKey, e);
            MessageUtils.showMessage("Error saving colors: " + e.getMessage());
            return false;
        }
    }

    /**
     * The highest priority scheme assigning colors to values of this attribute, or null.
     */
    private VariantColorScheme getSchemeForColors() {
        String key = infoKey.toUpperCase();
        for (VariantColorScheme scheme : VariantColorSchemes.getSchemes()) {
            if (scheme.covers(key)) {
                return scheme;
            }
        }
        return null;
    }

    private void repaintTrack() {
        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
    }
}
