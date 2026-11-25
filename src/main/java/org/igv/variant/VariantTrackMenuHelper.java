package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.sample.SampleMenuUtils;
import org.igv.sample.SampleSort;
import org.igv.track.AttributeManager;
import org.igv.track.DataType;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.renderer.ContinuousColorScale;
import org.igv.ui.IGV;
import org.igv.ui.util.MessageUtils;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
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


        items.addAll(TrackMenuUtils.getDisplayModeMenuItems(Collections.singletonList(variantTrack)));

        items.add(new JPopupMenu.Separator());
        items.add(new JLabel("<html>&nbsp;&nbsp;<b>Color variant sites by", JLabel.LEFT));
        items.add(getColorBandByAllelFrequency(variantTrack));
        items.add(getColorBandByAlleleFraction(variantTrack));
        JMenu infoFieldMenu = getColorByInfoFieldMenu(variantTrack);
        if (infoFieldMenu != null) {
            items.add(infoFieldMenu);
        }
        items.add(getEditColorItem(variantTrack));

        // Genotype color options
        if (variantTrack.isEnableMethylationRateSupport() || variantTrack.sampleCount() > 0) {
            items.add(new JPopupMenu.Separator());
            items.add(new JLabel("<html>&nbsp;&nbsp;<b>Color Samples By", JLabel.LEFT));
            items.add(getColorByGenotype(variantTrack));
            if (variantTrack.isEnableMethylationRateSupport()) {
                items.add(getColorByMethylationRate(variantTrack));
            }
            items.add(getColorByFormatFieldMenu(variantTrack));
            if (variantTrack.getColorGenotypesBy() != null) {
                items.add(getEditFormatColorItem(variantTrack));
            }
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
     * Opens the color dialog for what the track's variant sites are colored by: the INFO field, or the allele
     * frequency or fraction display -- a color scale or bars.
     */
    private static JMenuItem getEditColorItem(VariantTrack track) {

        final VariantTrack.ColorMode mode = track.getSiteColorMode();
        final String infoKey = track.getColorByAttribute();
        final boolean info = mode == VariantTrack.ColorMode.ATTRIBUTE && infoKey != null;

        JMenuItem item = new JMenuItem(info ? "INFO Field " + infoKey + " Colors..." : "Allele Frequency/Fraction Colors...");
        item.addActionListener(evt -> {
            Frame frame = IGV.getInstance().getMainFrame();
            if (info) {
                new VariantColorLegendDialog(frame, track, infoKey).setVisible(true);
            } else if (track.isAlleleFrequencyBars()) {
                new AlleleFrequencyBarColorsDialog(frame, track).setVisible(true);
            } else {
                String colorBy = mode == VariantTrack.ColorMode.ALLELE_FRACTION ?
                        VariantTrack.COLOR_BY_ALLELE_FRACTION : VariantTrack.COLOR_BY_ALLELE_FREQUENCY;
                new AlleleFrequencyColorsDialog(frame, track, colorBy).setVisible(true);
            }
        });
        return item;
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

        // A JMenu has no check box -- show one (checked when a field is selected, empty if the look-and-feel draws
        // one for unselected items), so the selection is visible without opening the submenu
        boolean selected = track.getSiteColorMode() == VariantTrack.ColorMode.ATTRIBUTE && track.getColorByAttribute() != null;
        JMenu menu = new CheckedMenu("INFO Field", selected);

        List<JMenuItem> fieldItems = new ArrayList<>(infoFields.size());
        for (VCFInfoHeaderLine infoField : infoFields) {
            fieldItems.add(getColorByInfoFieldItem(track, infoField));
        }

        addFieldItems(menu, fieldItems);

        return menu;
    }

    /**
     * Add the items for a file's INFO or FORMAT fields to a menu -- directly, or if there are too many for a
     * single menu, in alphabetical groups named for the first and last field of each.
     */
    private static void addFieldItems(JMenu menu, List<JMenuItem> fieldItems) {
        if (fieldItems.size() <= MAX_FLAT_INFO_FIELDS) {
            fieldItems.forEach(menu::add);
        } else {
            for (int i = 0; i < fieldItems.size(); i += INFO_FIELD_GROUP_SIZE) {
                int end = Math.min(i + INFO_FIELD_GROUP_SIZE, fieldItems.size());
                JMenu group = new JMenu(fieldItems.get(i).getText() + " - " + fieldItems.get(end - 1).getText());
                fieldItems.subList(i, end).forEach(group::add);
                menu.add(group);
            }
        }
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
     * The first time one is selected, give the track a default color scale over the loaded values, without asking.
     * Attributes the track or a scheme already covers, and attributes that are not numeric, need nothing.
     *
     * @return true if the track can be colored by this attribute
     */
    // Package private for testing
    static boolean defineScaleIfNeeded(VariantTrack track, String infoKey) {

        // Any scheme covering the attribute settles it -- with a scale, a categorical declaration, or discrete
        // colors.  A user who wrote discrete colors for a numeric attribute meant it.
        // A well-known allele frequency field has a scale already -- it is colored by rarity.
        // A track with a scale for the attribute, or colors for its values -- it was chosen on this track to be
        // categorical -- is settled too.
        if (!track.isNumericAttribute(infoKey) || VariantColorSchemes.getKeys().contains(infoKey.toUpperCase())
                || track.getAttributeScale(infoKey) != null || !track.getAllAttributeColors(infoKey).isEmpty()) {
            return true;
        }

        // Copy number has a scale of its own, whatever values happen to be loaded
        if (VariantTrack.formatKey("CN").equals(infoKey)) {
            track.setColorScale(infoKey, copyNumberScale());
            return true;
        }

        double[] range = track.getAttributeRange(infoKey);

        if (range == null) {
            // No numeric values loaded -- none, or none that can be read as numbers, such as multi-valued records.
            // A scale over a made-up range would color everything alike, and coloring by value anyway would quietly
            // treat a quantity as a category, so the same click would behave differently after navigating.  Leave
            // the selection alone and say why.
            String name = VariantTrack.attributeName(infoKey);
            MessageUtils.showMessage("No loaded features have a numeric value for " + name + ".  Move to a region where "
                    + name + " has values and choose it again"
                    // Color schemes are for INFO attributes
                    + (VariantTrack.isFormatKey(infoKey) ? "." : ", or import a color scheme for " + name + "."));
            return false;
        }

        ContinuousColorScale scale =
                new ContinuousColorScale(range[0], rangeEnd(range), DEFAULT_SCALE_MIN_COLOR, DEFAULT_SCALE_MAX_COLOR);

        track.setColorScale(infoKey, scale);
        return true;
    }

    /**
     * The default scale for the FORMAT field CN, the copy number genotype: the colors of the scale segmented copy
     * number tracks use.  Those hold log2(copy number / 2) -- neutral at 0 and saturated at +/- 1.5 -- which for a
     * copy number is neutral at 2 and saturated at 0 and at 6 (2 * 2^1.5 is 5.7).
     */
    // Package private for testing
    static ContinuousColorScale copyNumberScale() {
        ContinuousColorScale segScale = PreferencesManager.getPreferences().getColorScale(DataType.COPY_NUMBER);
        return new ContinuousColorScale(2, 0, 2, 6, segScale.getMinColor(), segScale.getMidColor(), segScale.getMaxColor());
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

        if (track.getAttributeScale(infoKey) != null) {
            return true;
        }

        int count = track.getAttributeValues(infoKey).size();
        if (count <= MAX_CATEGORICAL_VALUES) {
            return true;
        }

        MessageUtils.showMessage("Color by is not available for attributes with more than " + MAX_CATEGORICAL_VALUES
                + " distinct values.  " + VariantTrack.attributeName(infoKey) + " has " + count + " among the loaded features.");
        return false;
    }

    /**
     * A scale over a single repeated value has nothing to shade across, so give it somewhere to go.
     */
    private static double rangeEnd(double[] range) {
        return range[1] > range[0] ? range[1] : range[0] + 1;
    }

    private static JMenuItem getShowGenotypes(VariantTrack track) {
        // Reflects what is drawn -- a collapsed track shows no genotypes, and turning them on expands it
        final JMenuItem item = new JCheckBoxMenuItem("Show Genotypes", track.areGenotypesShown());
        item.addActionListener(evt -> {
            track.setGenotypesShown(item.isSelected());
            IGV.getInstance().revalidateTrackPanels();
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    /**
     * Menu for coloring genotypes by a VCF FORMAT field.  A field with no values among the loaded features is
     * disabled, as there is nothing to color by -- and so is the menu if that is all of them.
     */
    private static JMenu getColorByFormatFieldMenu(VariantTrack track) {

        List<VCFFormatHeaderLine> formatFields = track.getColorableFormatFields();

        boolean selected = track.getColorGenotypesBy() != null;
        JMenu menu = new CheckedMenu("FORMAT Field", selected);

        Set<String> withValues = track.getFormatFieldsWithValues(
                formatFields.stream().map(VCFFormatHeaderLine::getID).toList());
        menu.setEnabled(selected || !withValues.isEmpty());
        if (!menu.isEnabled()) {
            menu.setToolTipText("No loaded features have values for a FORMAT field other than GT");
        }

        List<JMenuItem> fieldItems = new ArrayList<>(formatFields.size());
        for (VCFFormatHeaderLine formatField : formatFields) {
            JMenuItem item = getColorByFormatFieldItem(track, formatField);
            if (!withValues.contains(formatField.getID()) && !item.isSelected()) {
                item.setEnabled(false);
                item.setToolTipText("No loaded features have a value for " + formatField.getID());
            }
            fieldItems.add(item);
        }

        addFieldItems(menu, fieldItems);

        return menu;
    }

    private static JMenuItem getColorByFormatFieldItem(VariantTrack track, VCFFormatHeaderLine formatField) {
        final String id = formatField.getID();
        final JMenuItem item = new JCheckBoxMenuItem(id, id.equals(track.getColorGenotypesBy()));
        String description = formatField.getDescription();
        if (description != null && description.length() > 0) {
            item.setToolTipText(description);
        }
        item.addActionListener(evt -> {
            String key = VariantTrack.formatKey(id);
            if (defineScaleIfNeeded(track, key) && isWithinCategoryLimit(track, key)) {
                track.setColorGenotypesBy(id);
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }

    /**
     * Opens the color dialog for the FORMAT field the track's genotypes are colored by.
     */
    private static JMenuItem getEditFormatColorItem(VariantTrack track) {
        final String id = track.getColorGenotypesBy();
        JMenuItem item = new JMenuItem("FORMAT Field " + id + " Colors...");
        item.addActionListener(evt ->
                new VariantColorLegendDialog(IGV.getInstance().getMainFrame(), track, VariantTrack.formatKey(id)).setVisible(true));
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


    /**
     * A menu drawn with the check box of a check box menu item -- checked or not, as the look-and-feel draws it --
     * so a selection shows without opening the menu.  A JMenu has no check box, and an icon would go in the icon
     * column, not the check column, under some look-and-feels (Aqua, Metal, Windows).  So the mark is taken from
     * the look-and-feel itself: a selected and an unselected check box item of this menu's size and text are
     * painted, and the pixels that differ are drawn over the menu.  The text is needed: some look-and-feels (Aqua)
     * position the mark by it.  Some look-and-feels (Metal) also draw a box that is the same selected or not; it
     * is found by comparing an unselected check box item with a plain item, both without text.
     */
    private static class CheckedMenu extends JMenu {

        private final boolean checked;

        CheckedMenu(String text, boolean checked) {
            super(text);
            this.checked = checked;
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            int w = getWidth();
            int h = getHeight();
            if (w <= 0 || h <= 0) {
                return;
            }
            BufferedImage selected = paintItem(new JCheckBoxMenuItem(getText(), true), w, h);
            BufferedImage unselected = paintItem(new JCheckBoxMenuItem(getText(), false), w, h);
            BufferedImage emptyCheck = paintItem(new JCheckBoxMenuItem("", false), w, h);
            BufferedImage emptyPlain = paintItem(new JMenuItem(""), w, h);
            BufferedImage mark = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    if (checked && selected.getRGB(x, y) != unselected.getRGB(x, y)) {
                        mark.setRGB(x, y, selected.getRGB(x, y));
                    } else if (emptyCheck.getRGB(x, y) != emptyPlain.getRGB(x, y)) {
                        mark.setRGB(x, y, emptyCheck.getRGB(x, y));
                    }
                }
            }
            g.drawImage(mark, 0, 0, null);
        }

        private BufferedImage paintItem(JMenuItem item, int w, int h) {
            item.setFont(getFont());
            item.setEnabled(isEnabled());
            // Highlighted when this menu is -- under the mouse, or open
            item.getModel().setArmed(getModel().isArmed() || isSelected());
            item.setSize(w, h);
            BufferedImage image = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
            Graphics2D g = image.createGraphics();
            item.paint(g);
            g.dispose();
            return image;
        }
    }
}
