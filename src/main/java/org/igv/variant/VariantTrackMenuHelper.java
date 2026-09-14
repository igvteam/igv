package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.sample.SampleMenuUtils;
import org.igv.sample.SampleSort;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.renderer.ContinuousColorScale;
import org.igv.ui.IGV;
import org.igv.ui.legend.HeatmapLegendEditor;
import org.igv.ui.util.MessageUtils;

import htsjdk.variant.vcf.VCFInfoHeaderLine;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * User: Jesse Whitworth
 * Date: Jul 16, 2010
 */
public class VariantTrackMenuHelper {

    private static Logger log = LogManager.getLogger(VariantTrackMenuHelper.class);
    private VariantTrack track;
    private static boolean depthSortingDirection;
    private static boolean genotypeSortingDirection;
    private static boolean sampleSortingDirection;
    private static boolean qualitySortingDirection;

    /**
     * Maximum number of INFO fields listed directly in the "color by" menu.  Beyond this they are broken into
     * alphabetical submenus -- some annotation pipelines define hundreds of them.
     */
    private static final int MAX_FLAT_INFO_FIELDS = 25;

    /**
     * Starting colors offered for a new color scale, pale to saturated.
     */
    private static final Color DEFAULT_SCALE_MIN_COLOR = new Color(255, 255, 204);
    private static final Color DEFAULT_SCALE_MAX_COLOR = new Color(202, 0, 32);

    private static final int INFO_FIELD_GROUP_SIZE = 20;

    /**
     * Most distinct values an attribute can have and still be colored by category.  Colors stop being
     * distinguishable, or useful, in the tens of values; past this the legend is unusable and the palette would
     * be handing out hundreds of colors.
     */
    static final int MAX_CATEGORICAL_VALUES = 500;

    /**
     * Return menu items for the variant track popup menu.
     */
    static List<Component> getMenuItems(final VariantTrack variantTrack, final Variant variant, TrackClickEvent e) {

        List<Component> items = new ArrayList<>();


        items.addAll(TrackMenuUtils.getSquishExpandItems(Collections.singletonList(variantTrack)));
        items.add(TrackMenuUtils.getRowHeightItem(Collections.singletonList(variantTrack)));
        items.add(TrackMenuUtils.getMinimizeHeightItem(Collections.singletonList(variantTrack)));

        items.add(new JPopupMenu.Separator());
        items.add(new JLabel("<html>&nbsp;&nbsp;<b>Color By", JLabel.LEFT));
        items.add(getColorBandByAllelFrequency(variantTrack));
        items.add(getColorBandByAlleleFraction(variantTrack));
        items.add(getAlleleFrequencyDisplayMenu(variantTrack));
        JMenu infoFieldMenu = getColorByInfoFieldMenu(variantTrack);
        if (infoFieldMenu != null) {
            items.add(infoFieldMenu);
            items.add(getEditInfoColorsItem(variantTrack));
        }
        items.add(getColorByNone(variantTrack));

        // Methylation color options
        if (variantTrack.isEnableMethylationRateSupport()) {
            items.add(new JPopupMenu.Separator());
            items.add(new JLabel("<html>&nbsp;&nbsp;<b>Color Samples By", JLabel.LEFT));
            items.add(getColorByGenotype(variantTrack));
            items.add(getColorByMethylationRate(variantTrack));
        }

        // Show genotypes
        if (variantTrack.sampleCount() > 0) {
            items.add(new JPopupMenu.Separator());
            items.add(getShowGenotypes(variantTrack));
        }

        //Sorter
        items.add(new JPopupMenu.Separator());
        for (JMenuItem item : getSortMenuItems(variantTrack, variant)) {
            items.add(item);
            item.setEnabled(variant != null);
        }

        boolean hasAttributes = AttributeManager.getInstance().getVisibleAttributes().size() > 0;
        if (hasAttributes || variantTrack.hasSamples()) {
            items.add(new JPopupMenu.Separator());
            if (hasAttributes) {
                items.add(SampleMenuUtils.getSortByAttributeItem(variantTrack));
                items.add(SampleMenuUtils.getGroupByAttributeItem(variantTrack));
                items.add(SampleMenuUtils.getFilterByAttributeItem(variantTrack));
            }
            if (variantTrack.hasSamples()) {
                items.add(SampleMenuUtils.getFilterByIdItem(variantTrack));
            }
        }

        items.add(new JPopupMenu.Separator());
        JMenuItem circItem = new JMenuItem("Add SVs to Circular View");
        circItem.addActionListener(e1 -> variantTrack.sendToCircularView(e));
        items.add(circItem);
        items.add(new JPopupMenu.Separator());

        items.add(getHideFilteredItem(variantTrack));
        items.add(getFeatureVisibilityItem(variantTrack));

        return items;
    }

    private static JMenuItem getFeatureVisibilityItem(VariantTrack track) {
        JMenuItem item = new JMenuItem("Set Feature Visibility Window...");
        item.addActionListener(evt -> {
            changeVisibilityWindow(track);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    private static JMenuItem getColorBandByAllelFrequency(VariantTrack track) {
        final JMenuItem item = new JCheckBoxMenuItem("Allele Frequency", track.getSiteColorMode() == VariantTrack.ColorMode.ALLELE_FREQUENCY);
        item.addActionListener(evt -> {
            track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    private static JMenuItem getColorBandByAlleleFraction(VariantTrack track) {
        final JMenuItem item = new JCheckBoxMenuItem("Allele Fraction", track.getSiteColorMode() == VariantTrack.ColorMode.ALLELE_FRACTION);
        item.addActionListener(evt -> {
            track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FRACTION);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    /**
     * How allele frequency and fraction are drawn: by a color scale, or as a bar whose height is the frequency.
     * Choosing a display while coloring by something else switches to allele frequency, so the choice is visible.
     */
    private static JMenu getAlleleFrequencyDisplayMenu(VariantTrack track) {

        JMenu menu = new JMenu("Allele Frequency Display");
        boolean bars = track.isAlleleFrequencyBars();

        JRadioButtonMenuItem colorItem = new JRadioButtonMenuItem("Color Scale", !bars);
        colorItem.addActionListener(evt -> setAlleleFrequencyDisplay(track, false));
        JRadioButtonMenuItem barItem = new JRadioButtonMenuItem("Bar Height", bars);
        barItem.addActionListener(evt -> setAlleleFrequencyDisplay(track, true));

        ButtonGroup group = new ButtonGroup();
        group.add(colorItem);
        group.add(barItem);
        menu.add(colorItem);
        menu.add(barItem);

        menu.addSeparator();
        // The colors, and legend, for the display the track shows
        JMenuItem editItem = new JMenuItem("Allele Frequency Colors...");
        editItem.addActionListener(evt -> {
            Frame frame = IGV.getInstance().getMainFrame();
            if (track.isAlleleFrequencyBars()) {
                new AlleleFrequencyBarColorsDialog(frame, track).setVisible(true);
            } else {
                new AlleleFrequencyColorsDialog(frame).setVisible(true);
            }
        });
        menu.add(editItem);

        return menu;
    }

    private static void setAlleleFrequencyDisplay(VariantTrack track, boolean bars) {
        track.setAlleleFrequencyBars(bars);
        VariantTrack.ColorMode mode = track.getSiteColorMode();
        if (mode != VariantTrack.ColorMode.ALLELE_FREQUENCY && mode != VariantTrack.ColorMode.ALLELE_FRACTION) {
            track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
        }
        IGV.getInstance().getContentPane().repaint();
    }

    /**
     * Menu for coloring the variant band by a VCF INFO attribute.  Returns null if the file has no INFO fields
     * that can be colored by.
     */
    private static JMenu getColorByInfoFieldMenu(VariantTrack track) {

        List<VCFInfoHeaderLine> infoFields = getMenuInfoFields(track);
        if (infoFields.isEmpty()) {
            return null;
        }

        JMenu menu = new JMenu("INFO Field");

        List<JMenuItem> fieldItems = new ArrayList<>(infoFields.size());
        for (VCFInfoHeaderLine infoField : infoFields) {
            fieldItems.add(getColorByInfoFieldItem(track, infoField));
        }

        if (fieldItems.size() <= MAX_FLAT_INFO_FIELDS) {
            fieldItems.forEach(menu::add);
        } else {
            // Too many fields for a single menu -- break them into alphabetical groups
            for (int i = 0; i < fieldItems.size(); i += INFO_FIELD_GROUP_SIZE) {
                int end = Math.min(i + INFO_FIELD_GROUP_SIZE, fieldItems.size());
                JMenu group = new JMenu(infoFields.get(i).getID() + " - " + infoFields.get(end - 1).getID());
                fieldItems.subList(i, end).forEach(group::add);
                menu.add(group);
            }
        }

        return menu;
    }

    /**
     * The INFO fields offered in the color-by menu.  AF is left out, as it has its own item, "Allele Frequency".
     */
    // Package private for testing
    static List<VCFInfoHeaderLine> getMenuInfoFields(VariantTrack track) {
        return track.getColorableInfoFields().stream()
                .filter(line -> !"AF".equals(line.getID()))
                .toList();
    }

    private static JMenuItem getColorByInfoFieldItem(VariantTrack track, VCFInfoHeaderLine infoField) {
        final String id = infoField.getID();
        final JMenuItem item = new JCheckBoxMenuItem(id,
                track.getSiteColorMode() == VariantTrack.ColorMode.ATTRIBUTE && id.equals(track.getColorByAttribute()));
        String description = infoField.getDescription();
        if (description != null && description.length() > 0) {
            item.setToolTipText(description);
        }
        item.addActionListener(evt -> {
            if (defineScaleIfNeeded(track, id) && isWithinCategoryLimit(track, id)) {
                track.setColorByAttribute(id);
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }

    /**
     * A numeric attribute is a quantity, not a category -- coloring it by value would give a color per variant.
     * The first time one is selected, have the user define a color scale, which is saved as a scheme so it
     * applies to every VCF with this attribute and can be shared.  Attributes a scheme already covers, and
     * attributes that are not numeric, need nothing.
     *
     * @return true if the track can be colored by this attribute
     */
    // Package private for testing -- the paths that do not open the dialog are worth covering
    static boolean defineScaleIfNeeded(VariantTrack track, String infoKey) {

        // Any scheme covering the attribute settles it -- with a scale, a categorical declaration, or discrete
        // colors.  A user who wrote discrete colors for a numeric attribute meant it.
        // A well-known allele frequency field has a scale already -- it is colored by rarity.
        if (!track.isNumericAttribute(infoKey) || VariantColorSchemes.getKeys().contains(infoKey.toUpperCase())
                || VariantColorSchemes.getScale(infoKey) != null) {
            return true;
        }

        double[] range = track.getAttributeRange(infoKey);

        if (range == null && track.getAttributeValues(infoKey).isEmpty()) {
            // Nothing loaded to judge by.  Coloring by value anyway would quietly treat a quantity as a category
            // and record nothing, so the same click would behave differently after navigating.  Leave the
            // selection alone and say why.
            MessageUtils.showMessage("No loaded features have a value for " + infoKey + ".  Move to a region where "
                    + infoKey + " has values and choose it again, or import a color scheme for " + infoKey + ".");
            return false;
        }

        // Values are loaded.  Prefill the range from them when they are numeric; when they are not -- every
        // record multi-valued, say -- still ask, rather than deciding the attribute is categorical unasked.
        ContinuousColorScale scale = range == null
                ? new ContinuousColorScale(0, 1, DEFAULT_SCALE_MIN_COLOR, DEFAULT_SCALE_MAX_COLOR)
                : new ContinuousColorScale(range[0], rangeEnd(range), DEFAULT_SCALE_MIN_COLOR, DEFAULT_SCALE_MAX_COLOR);

        // IGV cannot tell a quantity from a code, so the dialog offers both: define a scale, or say the values
        // are categories after all and colour them individually.
        HeatmapLegendEditor editor = new HeatmapLegendEditor(
                IGV.getInstance().getMainFrame(), true, scale, "Values Are Categories");
        editor.setTitle("Color scale for " + infoKey);
        editor.setVisible(true);

        if (editor.isDiscreteSelected()) {
            return saveDiscreteScheme(track, infoKey);
        }
        if (editor.isCanceled()) {
            return false;
        }

        try {
            VariantColorSchemes.saveScale(infoKey + " colorscale", infoKey, editor.getColorScheme());
            return true;
        } catch (Exception e) {
            log.error("Error saving color scale for " + infoKey, e);
            MessageUtils.showMessage("Error saving color scale: " + e.getMessage());
            return false;
        }
    }

    /**
     * An attribute colored by category may have at most MAX_CATEGORICAL_VALUES distinct values among the loaded
     * features; if it has more, say so and leave the track's coloring alone.  An attribute colored by a scale has
     * no limit -- a depth or a frequency can take any number of values.
     *
     * @return true if the track can be colored by this attribute
     */
    // Package private for testing
    static boolean isWithinCategoryLimit(VariantTrack track, String infoKey) {

        if (VariantColorSchemes.getScale(infoKey) != null) {
            return true;
        }

        int count = track.getAttributeValues(infoKey).size();
        if (count <= MAX_CATEGORICAL_VALUES) {
            return true;
        }

        MessageUtils.showMessage("Color by is not available for attributes with more than " + MAX_CATEGORICAL_VALUES
                + " distinct values.  " + infoKey + " has " + count + " among the loaded features.");
        return false;
    }

    /**
     * A scale over a single repeated value has nothing to shade across, so give it somewhere to go.
     */
    private static double rangeEnd(double[] range) {
        return range[1] > range[0] ? range[1] : range[0] + 1;
    }

    /**
     * Record that a numeric attribute holds categories, by saving a scheme with a color per loaded value.  The
     * scheme is what stops the scale dialog reappearing, and the colors are then editable from the legend and
     * shareable like any other.  Values seen later still get colors from the palette.
     */
    private static boolean saveDiscreteScheme(VariantTrack track, String infoKey) {

        // Checked before writing anything, so an attribute over the limit does not leave a scheme behind
        if (!isWithinCategoryLimit(track, infoKey)) {
            return false;
        }

        Map<String, Color> colors = new LinkedHashMap<>();
        for (String value : track.getAttributeValues(infoKey)) {
            colors.put(value, track.getAttributeColor(infoKey, value));
        }

        try {
            // The declaration persists the answer even when no values are loaded to color
            VariantColorSchemes.saveScheme(infoKey + " colors", infoKey, colors, true);
            return true;
        } catch (Exception e) {
            log.error("Error saving colors for " + infoKey, e);
            MessageUtils.showMessage("Error saving colors: " + e.getMessage());
            return false;
        }
    }

    /**
     * Opens the color legend for the attribute the track is colored by.  Only enabled in that mode -- there is
     * nothing to show a legend for otherwise.
     */
    private static JMenuItem getEditInfoColorsItem(VariantTrack track) {

        final String infoKey = track.getColorByAttribute();
        final boolean active = track.getSiteColorMode() == VariantTrack.ColorMode.ATTRIBUTE && infoKey != null;

        JMenuItem item = new JMenuItem(active ? "Edit " + infoKey + " Colors..." : "Edit INFO Colors...");
        item.setEnabled(active);
        if (active) {
            item.addActionListener(evt ->
                    new VariantColorLegendDialog(IGV.getInstance().getMainFrame(), track, infoKey).setVisible(true));
        }
        return item;
    }

    private static JMenuItem getColorByNone(VariantTrack track) {
        final JMenuItem item = new JCheckBoxMenuItem("None", track.getSiteColorMode() == VariantTrack.ColorMode.NONE);
        item.addActionListener(evt -> {
            track.setSiteColorMode(VariantTrack.ColorMode.NONE);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private static JMenuItem getShowGenotypes(VariantTrack track) {
        final JMenuItem item = new JCheckBoxMenuItem("Show Genotypes", track.isShowGenotypes());
        item.addActionListener(evt -> {
            track.setShowGenotypes(item.isSelected());
            IGV.getInstance().revalidateTrackPanels();
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private static JMenuItem getColorByGenotype(VariantTrack track) {
        final JMenuItem item = new JCheckBoxMenuItem("Genotype", track.getGenotypeColorMode() == VariantTrack.ColorMode.GENOTYPE);
        item.addActionListener(evt -> {
            track.setGenotypeColorMode(VariantTrack.ColorMode.GENOTYPE);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private static JMenuItem getColorByMethylationRate(VariantTrack track) {
        final JMenuItem item = new JCheckBoxMenuItem("Methylation Rate", track.getGenotypeColorMode() == VariantTrack.ColorMode.METHYLATION_RATE);
        item.addActionListener(evt -> {
            track.setGenotypeColorMode(VariantTrack.ColorMode.METHYLATION_RATE);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    private static JMenuItem getHideFilteredItem(VariantTrack track) {
        JMenuItem item = new JCheckBoxMenuItem("Hide Filtered Sites", track.getHideFiltered());
        item.addActionListener(evt -> {
            track.setHideFiltered(!track.getHideFiltered());
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    public static JMenuItem getGenotypeSortItem(VariantTrack track, final Variant variant) {

        JMenuItem item = new JMenuItem("Sort By Genotype");
        if (variant != null) {
            item.addActionListener(evt -> {
                track.sortSamples(SampleSort.GENOTYPE, variant, !genotypeSortingDirection);
                genotypeSortingDirection = !genotypeSortingDirection;
                IGV.getInstance().getContentPane().repaint();
            });
        }

        return item;
    }

    public static JMenuItem getSampleNameSortItem(VariantTrack track, final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Sample Name");
        if (variant != null) {
            item.addActionListener(evt -> {
                track.sortSamplesByName(sampleSortingDirection);
                sampleSortingDirection = !sampleSortingDirection;
                IGV.getInstance().getContentPane().repaint();
            });
        }
        return item;
    }

    public static JMenuItem getDepthSortItem(VariantTrack track, final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Depth");
        if (variant != null) {
            item.addActionListener(evt -> {
                track.sortSamples(SampleSort.DEPTH, variant, !depthSortingDirection);
                depthSortingDirection = !depthSortingDirection;
                IGV.getInstance().getContentPane().repaint();
            });

        }
        return item;
    }

    public static JMenuItem getQualitySortItem(VariantTrack track, final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Quality");
        if (variant != null) {
            double quality = variant.getPhredScaledQual();
            if (quality > -1) {
                item.addActionListener(evt -> {
                    track.sortSamples(SampleSort.QUALITY, variant, !qualitySortingDirection);
                    qualitySortingDirection = !qualitySortingDirection;
                    IGV.getInstance().getContentPane().repaint();
                });
            } else {
                item.setEnabled(false);
            }
        }

        return item;
    }

    public static void changeVisibilityWindow(VariantTrack track) {
        TrackMenuUtils.changeFeatureVisibilityWindow(Arrays.asList((Track) track));
    }

    public static Collection<JMenuItem> getSortMenuItems(VariantTrack track, Variant variant) {

        java.util.List<JMenuItem> items = new ArrayList<JMenuItem>();
        items.add(getGenotypeSortItem(track, variant));
        items.add(getSampleNameSortItem(track, variant));
        items.add(getDepthSortItem(track, variant));
        items.add(getQualitySortItem(track, variant));
        return items;
    }

}
