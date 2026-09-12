package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.sample.SampleMenuUtils;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.ui.IGV;

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
    private static final int INFO_FIELD_GROUP_SIZE = 20;

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

        if (AttributeManager.getInstance().getVisibleAttributes().size() > 0) {
            items.add(new JPopupMenu.Separator());
            items.add(SampleMenuUtils.getSortByAttributeItem(variantTrack));
            items.add(SampleMenuUtils.getGroupByAttributeItem(variantTrack));
            items.add(SampleMenuUtils.getFilterByAttributeItem(variantTrack));
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
     * Menu for coloring the variant band by a VCF INFO attribute.  Returns null if the file has no INFO fields
     * that can be colored by.
     */
    private static JMenu getColorByInfoFieldMenu(VariantTrack track) {

        List<VCFInfoHeaderLine> infoFields = track.getColorableInfoFields();
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

    private static JMenuItem getColorByInfoFieldItem(VariantTrack track, VCFInfoHeaderLine infoField) {
        final String id = infoField.getID();
        final JMenuItem item = new JCheckBoxMenuItem(id,
                track.getSiteColorMode() == VariantTrack.ColorMode.ATTRIBUTE && id.equals(track.getColorByAttribute()));
        String description = infoField.getDescription();
        if (description != null && description.length() > 0) {
            item.setToolTipText(description);
        }
        item.addActionListener(evt -> {
            track.setColorByAttribute(id);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
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
                track.sortSamples(new GenotypeComparator(variant, genotypeSortingDirection));
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
                Comparator<String> comparator = sampleSortingDirection ? String::compareTo : (s1, s2) -> s2.compareTo(s1);
                track.sortSamples(comparator);
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
                track.sortSamples(new DepthComparator(variant, depthSortingDirection));
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
                    track.sortSamples(new QualityComparator(variant, qualitySortingDirection));
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
