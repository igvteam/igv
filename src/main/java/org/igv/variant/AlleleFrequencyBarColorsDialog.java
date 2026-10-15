package org.igv.variant;

import org.igv.Globals;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.color.ColorSwatch;
import org.igv.ui.color.ColorUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

/**
 * Legend and colors for the allele frequency bar display.  The lower part of a variant's bar, in the variant color, is
 * the frequency of the alternate alleles; the rest is the reference color.  The colors default to the AF_VAR.COLOR
 * and AF_REF.COLOR preferences (Preferences > Variants).  Changes apply to the track, and with "Update Preferences"
 * checked are also saved as those preferences.
 */
public class AlleleFrequencyBarColorsDialog extends JDialog {

    /**
     * As drawn by {@link VariantRenderer} for a reference block (only a &lt;NON_REF&gt; or &lt;*&gt; alternate allele)
     */
    private static final Color REFERENCE_BLOCK_COLOR = new Color(200, 200, 215);

    private final VariantTrack track;
    private final IGVPreferences prefs = PreferencesManager.getPreferences();

    private final JCheckBox updatePreferencesBox = new JCheckBox("Update Preferences");

    private Color variantColor;
    private Color referenceColor;
    private boolean variantColorChanged;
    private boolean referenceColorChanged;
    private boolean resetToPreferences;

    private final ColorSwatch variantSwatch;
    private final ColorSwatch referenceSwatch;
    private final SampleBars sampleBars = new SampleBars();

    public AlleleFrequencyBarColorsDialog(Frame owner, VariantTrack track) {
        super(owner, "Allele Frequency Colors", true);
        this.track = track;

        variantColor = track.getAlleleFrequencyVariantColor();
        referenceColor = track.getAlleleFrequencyReferenceColor();

        variantSwatch = new ColorSwatch(variantColor);
        variantSwatch.addColorChangeListener(c -> {
            variantColor = c;
            variantColorChanged = true;
            sampleBars.repaint();
        });
        referenceSwatch = new ColorSwatch(referenceColor);
        referenceSwatch.addColorChangeListener(c -> {
            referenceColor = c;
            referenceColorChanged = true;
            sampleBars.repaint();
        });

        initComponents();
        pack();
        setLocationRelativeTo(owner);
    }

    private void initComponents() {

        setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                cancel();
            }
        });

        JPanel content = new JPanel(new BorderLayout(0, 8));
        content.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        setContentPane(content);

        String value = track.getSiteColorMode() == VariantTrack.ColorMode.ALLELE_FRACTION ?
                "the allele fraction: AC divided by AN" :
                "the allele frequency: AF, or GMAF if a variant has no AF";
        content.add(new JLabel("<html>The lower part of each bar, in the variant color, is " + value + ".<br>" +
                "At a site with more than one alternate allele, the alleles are added.<br>" +
                "A variant with no value, or a value of 0, is drawn in the track color, " +
                "and a gVCF reference block in gray.<br>" +
                "Filtered variants are drawn faded.  Changes apply to this track."), BorderLayout.NORTH);

        JPanel swatches = new JPanel();
        swatches.setLayout(new BoxLayout(swatches, BoxLayout.Y_AXIS));
        swatches.add(swatchRow(variantSwatch, "Variant (alternate alleles)"));
        swatches.add(swatchRow(referenceSwatch, "Reference"));

        JPanel center = new JPanel(new BorderLayout(0, 8));
        center.add(sampleBars, BorderLayout.NORTH);
        center.add(swatches, BorderLayout.CENTER);
        content.add(center, BorderLayout.CENTER);

        JButton resetButton = new JButton("Reset to Preferences");
        resetButton.setToolTipText("Show the colors from preferences -- with \"OK\", the track follows them again");
        resetButton.addActionListener(e -> resetToPreferences());
        // There is always a preference -- the default colors if the user never set them -- so it is not replaced
        // unless asked
        updatePreferencesBox.setToolTipText("Also save these colors as the preference for all variant tracks");
        JPanel leftButtons = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        leftButtons.add(resetButton);
        leftButtons.add(updatePreferencesBox);

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> save());
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> cancel());
        JPanel rightButtons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 4, 0));
        rightButtons.add(okButton);
        rightButtons.add(cancelButton);

        JPanel buttons = new JPanel(new BorderLayout());
        buttons.add(leftButtons, BorderLayout.WEST);
        buttons.add(rightButtons, BorderLayout.EAST);
        content.add(buttons, BorderLayout.SOUTH);

        getRootPane().setDefaultButton(okButton);
    }

    private static JPanel swatchRow(ColorSwatch swatch, String label) {
        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));
        row.add(swatch);
        row.add(new JLabel(label));
        return row;
    }

    /**
     * Show the preference colors.  With "OK", and no edit since, the track follows the preferences again.
     */
    private void resetToPreferences() {
        variantColor = prefs.getAsColor(Constants.AF_VAR_COLOR);
        referenceColor = prefs.getAsColor(Constants.AF_REF_COLOR, Globals.DARK_MODE_BLUE);
        variantColorChanged = false;
        referenceColorChanged = false;
        resetToPreferences = true;
        variantSwatch.setSelectedColor(variantColor);
        referenceSwatch.setSelectedColor(referenceColor);
        sampleBars.repaint();
    }

    private void save() {

        if (updatePreferencesBox.isSelected()) {
            if (variantColorChanged) {
                prefs.put(Constants.AF_VAR_COLOR, ColorUtilities.colorToString(variantColor));
            }
            if (referenceColorChanged) {
                prefs.put(Constants.AF_REF_COLOR, ColorUtilities.colorToString(referenceColor));
            }
        }

        // The track keeps a color it was given, edited here or before, unless reset to the preference
        if (resetToPreferences) {
            track.setAlleleFrequencyVariantColor(null);
            track.setAlleleFrequencyReferenceColor(null);
        }
        if (variantColorChanged) {
            track.setAlleleFrequencyVariantColor(variantColor);
        }
        if (referenceColorChanged) {
            track.setAlleleFrequencyReferenceColor(referenceColor);
        }

        repaintTracks();
        dispose();
    }

    private void cancel() {
        dispose();
    }

    private static void repaintTracks() {
        if (IGV.hasInstance()) {
            IGV.getInstance().repaint();
        }
    }

    /**
     * Example bars: two frequencies, a variant with no value (or 0), and a reference block.
     */
    private class SampleBars extends JComponent {

        private static final int BAR_WIDTH = 24;
        private static final int BAR_HEIGHT = 40;
        private static final int SPACING = 100;

        SampleBars() {
            setPreferredSize(new Dimension(4 * SPACING, BAR_HEIGHT + 24));
        }

        @Override
        protected void paintComponent(Graphics g) {

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setFont(g2.getFont().deriveFont(g2.getFont().getSize() - 1f));
            FontMetrics fm = g2.getFontMetrics();

            drawBar(g2, fm, 0, 0.3, "30%", null);
            drawBar(g2, fm, 1, 0.8, "80%", null);
            drawBar(g2, fm, 2, 0, "No value or 0%", track.getColor());
            drawBar(g2, fm, 3, 0, "Reference block", REFERENCE_BLOCK_COLOR);
            g2.dispose();
        }

        /**
         * @param fill a solid color for the whole bar, or null to split it by frequency
         */
        private void drawBar(Graphics2D g2, FontMetrics fm, int index, double frequency, String label, Color fill) {

            int center = SPACING / 2 + index * SPACING;
            int x = center - BAR_WIDTH / 2;

            if (fill != null) {
                g2.setColor(fill);
                g2.fillRect(x, 0, BAR_WIDTH, BAR_HEIGHT);
            } else {
                int variantHeight = (int) Math.round(frequency * BAR_HEIGHT);
                g2.setColor(referenceColor);
                g2.fillRect(x, 0, BAR_WIDTH, BAR_HEIGHT - variantHeight);
                g2.setColor(variantColor);
                g2.fillRect(x, BAR_HEIGHT - variantHeight, BAR_WIDTH, variantHeight);
            }

            g2.setColor(getForeground());
            g2.drawString(label, center - fm.stringWidth(label) / 2, BAR_HEIGHT + 4 + fm.getAscent());
        }
    }
}
