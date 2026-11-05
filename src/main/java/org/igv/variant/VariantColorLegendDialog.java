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
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
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
 * at which the user can see what needs recoloring.  Nothing changes until "OK", which applies the edits to this
 * track -- and, with "Update Preferences" checked, also saves them in a scheme that applies to every VCF with
 * this attribute.  "Cancel" discards them.
 */
public class VariantColorLegendDialog extends JDialog {

    private static Logger log = LogManager.getLogger(VariantColorLegendDialog.class);

    private final VariantTrack track;
    private final String infoKey;
    private final JPanel valuePanel = new JPanel();
    private final JLabel description = new JLabel();
    private final JButton editScaleButton = new JButton("Edit Scale...");
    private final JButton categoriesButton = new JButton("Values Are Categories...");
    private final JCheckBox updatePreferencesBox = new JCheckBox("Update Preferences");

    // Edits not yet saved -- applied on "OK", discarded on "Cancel"
    private final Map<String, Color> pendingColors = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);
    private ContinuousColorScale pendingScale;

    // The track's color table, read once per populate() -- it can hold thousands of values
    private Map<String, Color> trackColors = Map.of();

    public VariantColorLegendDialog(Frame owner, VariantTrack track, String infoKey) {

        super(owner, "INFO Field " + infoKey + " Colors", true);
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
        categoriesButton.setToolTipText("Color each value of this numeric attribute separately, instead of by a scale");
        categoriesButton.addActionListener(e -> useCategories());

        // Checked when there is no preference for the attribute yet, so the first colors chosen become the default.
        // An existing preference is not replaced unless asked.
        updatePreferencesBox.setSelected(VariantColorSchemes.getScheme(infoKey) == null);
        updatePreferencesBox.setToolTipText("Also save these colors so they apply to every VCF with this attribute");

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> {
            if (save()) {
                dispose();
            }
        });

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> dispose());

        JPanel leftButtons = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        leftButtons.add(updatePreferencesBox);
        JPanel rightButtons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 0));
        rightButtons.add(editScaleButton);
        rightButtons.add(categoriesButton);
        rightButtons.add(okButton);
        rightButtons.add(cancelButton);
        JPanel buttonPanel = new JPanel(new BorderLayout());
        buttonPanel.add(leftButtons, BorderLayout.WEST);
        buttonPanel.add(rightButtons, BorderLayout.EAST);
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
        trackColors = track.getAttributeColors(infoKey);

        AbstractColorScale scale = pendingScale != null ? pendingScale : track.getAttributeScale(infoKey);

        description.setText("<html>Values of <b>" + infoKey + "</b> among the loaded features.  "
                + (scale == null ? "Click a color to change it." : "Colored by a scale."));

        editScaleButton.setVisible(scale instanceof ContinuousColorScale || scale instanceof ColorStopScale);
        // A numeric attribute colored by a scale can instead be said to hold categories
        categoriesButton.setVisible(scale instanceof ContinuousColorScale && track.isNumericAttribute(infoKey));
        // Allele frequency colors are edited, and saved, in their own dialog
        updatePreferencesBox.setVisible(!(scale instanceof ColorStopScale));

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
     * The values to list: those among the loaded features (see {@link VariantTrack#getAttributeValues}), plus the
     * track's "*" entry if it has one.  Sorted case insensitively, keeping the spelling seen in the loaded features.
     * <p>
     * Values in the track's color table that are no longer loaded are deliberately left out.  The table grows with
     * every value ever drawn -- thousands, for an attribute like a gene name.  Their colors are not lost: the track
     * keeps them, they are saved in sessions, and they come back when the values do.
     */
    private Set<String> valuesToShow() {
        Set<String> values = new TreeSet<>(String.CASE_INSENSITIVE_ORDER);
        values.addAll(track.getAttributeValues(infoKey));
        if (trackColors.containsKey(VariantTrack.WILDCARD_VALUE)) {
            values.add(VariantTrack.WILDCARD_VALUE);
        }
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
        if (VariantTrack.WILDCARD_VALUE.equals(value)) {
            return "(all other values)";
        }
        if (trackColors.containsKey(value)) {
            return "(this track)";
        }
        VariantColorScheme scheme = VariantColorSchemes.getScheme(infoKey);
        if (scheme != null && scheme.getColor(infoKey, value) != null) {
            return "(" + scheme.getName() + ")";
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

        AbstractColorScale scale = pendingScale != null ? pendingScale : track.getAttributeScale(infoKey);
        if (scale instanceof ColorStopScale) {
            // A well-known allele frequency field -- edited, and saved, like the allele frequency modes
            new AlleleFrequencyColorsDialog(getOwner() instanceof Frame frame ? frame : null, track, infoKey)
                    .setVisible(true);
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
     * Say that a numeric attribute holds categories: drop its scale and give each loaded value its own color.  Applied
     * at once, not on "OK".  A scheme with a scale for the attribute would still win, so it is replaced, as it is when
     * "Update Preferences" is checked.
     */
    private void useCategories() {

        List<String> values = new ArrayList<>(track.getAttributeValues(infoKey));
        if (values.size() > VariantTrackMenuHelper.MAX_CATEGORICAL_VALUES) {
            MessageUtils.showMessage("Color by is not available for attributes with more than "
                    + VariantTrackMenuHelper.MAX_CATEGORICAL_VALUES + " distinct values.  " + infoKey + " has "
                    + values.size() + " among the loaded features.");
            return;
        }
        if (JOptionPane.showConfirmDialog(this, "Color each value of " + infoKey + " separately, instead of by a scale?",
                "Values Are Categories", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
            return;
        }

        List<Color> used = new ArrayList<>();
        Map<String, Color> colors = new LinkedHashMap<>();
        for (String value : values) {
            Color color = DistinctColors.next(used);
            used.add(color);
            colors.put(value, color);
        }

        if (updatePreferencesBox.isSelected() || VariantColorSchemes.getSchemeForScale(infoKey) != null) {
            try {
                // The declaration persists the answer even when no values are loaded to color
                VariantColorSchemes.saveScheme(infoKey, colors, true);
            } catch (Exception e) {
                log.error("Error saving colors for " + infoKey, e);
                MessageUtils.showMessage("Error saving colors: " + e.getMessage());
                return;
            }
        }

        track.setColorScale(infoKey, null);
        colors.forEach((value, color) -> track.setAttributeColor(infoKey, value, color));
        pendingScale = null;
        pendingColors.clear();
        populate();
        pack();
        repaintTrack();
    }

    /**
     * Apply the edits to the track.  The track keeps them whether or not they are also saved as preferences, so a
     * session shows them to whoever opens it.
     *
     * @return true if applied, or there was nothing to apply
     */
    private boolean save() {

        if (pendingScale == null && pendingColors.isEmpty()) {
            return true;
        }
        if (updatePreferencesBox.isSelected() && !saveToPreferences()) {
            return false;
        }

        if (pendingScale != null) {
            track.setColorScale(infoKey, pendingScale);
        }
        pendingColors.forEach((value, color) -> track.setAttributeColor(infoKey, value, color));
        repaintTrack();
        return true;
    }

    /**
     * Save the edits into the scheme for the attribute, adding to or replacing its entries, so the change applies
     * to every track with the attribute.  A scheme shipped with IGV cannot be written, so a copy is saved, which
     * then replaces it.  Updating the existing scheme rather than writing a new one keeps what it holds beyond
     * these edits -- colors for other values, a categorical declaration.
     *
     * @return true if saved
     */
    private boolean saveToPreferences() {

        VariantColorScheme existing = pendingScale != null
                ? VariantColorSchemes.getSchemeForScale(infoKey) : VariantColorSchemes.getScheme(infoKey);
        VariantColorScheme updated = existing == null ? new VariantColorScheme(infoKey) : existing.copy();

        if (pendingScale != null) {
            updated.setScale(infoKey, pendingScale);
        } else {
            pendingColors.forEach((value, color) -> updated.setColor(infoKey, value, color));
            if (track.isNumericAttribute(infoKey)) {
                // Colored by value, so it holds categories -- without the declaration it would be asked for a scale
                updated.setCategorical(infoKey);
            }
        }

        try {
            VariantColorSchemes.save(updated);
            return true;
        } catch (Exception e) {
            log.error("Error saving colors for " + infoKey, e);
            MessageUtils.showMessage("Error saving colors: " + e.getMessage());
            return false;
        }
    }

    private void repaintTrack() {
        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
    }
}
