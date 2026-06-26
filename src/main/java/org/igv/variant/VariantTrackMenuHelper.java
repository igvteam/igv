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
     * Return menu items for the variant track popup menu.
     */
    static List<Component> getMenuItems(final VariantTrack variantTrack, final Variant variant, TrackClickEvent e) {

        List<Component> items = new ArrayList<>();


        items.add(TrackMenuUtils.getRowHeightItem(Collections.singletonList(variantTrack)));
        items.add(TrackMenuUtils.getMinimizeHeightItem(Collections.singletonList(variantTrack)));

        items.add(new JPopupMenu.Separator());
        items.add(new JLabel("<html>&nbsp;&nbsp;<b>Color By", JLabel.LEFT));
        items.add(getColorBandByAllelFrequency(variantTrack));
        items.add(getColorBandByAlleleFraction(variantTrack));

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
