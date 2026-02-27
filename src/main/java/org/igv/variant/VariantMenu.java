package org.igv.variant;

import org.igv.jbrowse.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.sample.AttributeComparator;
import org.igv.sample.SampleMenuUtils;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.ui.IGV;
import org.igv.ui.action.GroupSamplesMenuAction;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.util.SortDialog;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

/**
 * User: Jesse Whitworth
 * Date: Jul 16, 2010
 */
public class VariantMenu extends IGVPopupMenu {

    private static Logger log = LogManager.getLogger(VariantMenu.class);
    private VariantTrack track;
    private static boolean depthSortingDirection;
    private static boolean genotypeSortingDirection;
    private static boolean sampleSortingDirection;
    private static boolean qualitySortingDirection;

    public VariantMenu(final VariantTrack variantTrack, final Variant variant, TrackClickEvent e) {
        super();
        this.track = variantTrack;

        //Title
        JLabel popupTitle = new JLabel("<html><b>" + this.track.getName(), JLabel.LEFT);
        Font newFont = getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        add(popupTitle);


        if (PreferencesManager.getPreferences().getAsBoolean(Constants.CIRC_VIEW_ENABLED) && CircularViewUtilities.ping()) {
            addSeparator();
            JMenuItem circItem = new JMenuItem("Add SVs to Circular View");
            circItem.addActionListener(e1 -> track.sendToCircularView(e));
            add(circItem);
        }


        //Change Track Settings
        addSeparator();

        List<Track> selectedTracks = Arrays.asList(variantTrack);
        add(TrackMenuUtils.getTrackRenameItem(selectedTracks));
        add(TrackMenuUtils.getChangeFontSizeItem(selectedTracks));

        // Color items
//        addSeparator();
//        JMenuItem colorItem = new JMenuItem("Set Track Color...");
//        colorItem.addActionListener(evt -> TrackMenuUtils.changeTrackColor(selectedTracks));
//        add(colorItem);

        addSeparator();
        JLabel colorSiteByItem = new JLabel("<html>&nbsp;&nbsp;<b>Color By", JLabel.LEFT);
        add(colorSiteByItem);
        add(getColorBandByAllelFrequency());
        add(getColorBandByAlleleFraction());
        //add(getColorByNone());

        //Hides
        if (track.isEnableMethylationRateSupport()) {
            addSeparator();
            JLabel colorByItem = new JLabel("<html>&nbsp;&nbsp;<b>Color Samples By", JLabel.LEFT);
            add(colorByItem);
            add(getColorByGenotype());
            add(getColorByMethylationRate());
        }

        // Show genotypes
        if (track.sampleCount() > 0) {
            addSeparator();
            add(getShowGenotypes());
        }

        //Sorter
        addSeparator();
        for (JMenuItem item : getSortMenuItems(variant)) {
            add(item);
            item.setEnabled(variant != null);
        }

        if (AttributeManager.getInstance().getVisibleAttributes().size() > 0) {
            add(SampleMenuUtils.getSortByAttributeItem(track));
            addSeparator();
            add(getGenotypeGroupItem());
            add(SampleMenuUtils.getGroupByAttributeItem(track));
        }

        //Variant Information
        addSeparator();
        JLabel displayHeading = new JLabel("Display Mode", JLabel.LEFT);
        add(displayHeading);
        for (JMenuItem item : getDisplayModeItems()) {
            add(item);
        }

        addSeparator();
        JMenuItem item = new JMenuItem("Change Squished Row Height...");
        item.addActionListener(evt -> {
            int currentValue = track.getSquishedHeight();
            Integer newValue = TrackMenuUtils.getIntegerInput("Squished row height", currentValue);
            if (newValue != null) {
                track.setSquishedHeight(newValue);
                IGV.getInstance().getContentPane().repaint();
            }
        });
        add(item);

        add(getHideFilteredItem());
        add(getFeatureVisibilityItem());
    }

    private JMenuItem getFeatureVisibilityItem() {
        JMenuItem item = new JMenuItem("Set Feature Visibility Window...");
        item.addActionListener(evt -> {
            changeVisibilityWindow();
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    private JMenuItem getColorBandByAllelFrequency() {
        final JMenuItem item = new JCheckBoxMenuItem("Allele Frequency", track.getSiteColorMode() == VariantTrack.ColorMode.ALLELE_FREQUENCY);
        item.addActionListener(evt -> {
            track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FREQUENCY);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    private JMenuItem getColorBandByAlleleFraction() {
        final JMenuItem item = new JCheckBoxMenuItem("Allele Fraction", track.getSiteColorMode() == VariantTrack.ColorMode.ALLELE_FRACTION);
        item.addActionListener(evt -> {
            track.setSiteColorMode(VariantTrack.ColorMode.ALLELE_FRACTION);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private JMenuItem getColorByNone() {
        final JMenuItem item = new JCheckBoxMenuItem("None", track.getSiteColorMode() == VariantTrack.ColorMode.NONE);
        item.addActionListener(evt -> {
            track.setSiteColorMode(VariantTrack.ColorMode.NONE);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private JMenuItem getShowGenotypes() {
        final JMenuItem item = new JCheckBoxMenuItem("Show Genotypes", track.isShowGenotypes());
        item.addActionListener(evt -> {
            track.setShowGenotypes(item.isSelected());
            IGV.getInstance().revalidateTrackPanels();
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private JMenuItem getColorByGenotype() {
        final JMenuItem item = new JCheckBoxMenuItem("Genotype", track.getGenotypeColorMode() == VariantTrack.ColorMode.GENOTYPE);
        item.addActionListener(evt -> {
            track.setGenotypeColorMode(VariantTrack.ColorMode.GENOTYPE);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }

    private JMenuItem getColorByMethylationRate() {
        final JMenuItem item = new JCheckBoxMenuItem("Methylation Rate", track.getGenotypeColorMode() == VariantTrack.ColorMode.METHYLATION_RATE);
        item.addActionListener(evt -> {
            track.setGenotypeColorMode(VariantTrack.ColorMode.METHYLATION_RATE);
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    private JMenuItem getHideFilteredItem() {
        JMenuItem item = new JCheckBoxMenuItem("Hide Filtered Sites", track.getHideFiltered());
        item.addActionListener(evt -> {
            track.setHideFiltered(!track.getHideFiltered());
            IGV.getInstance().getContentPane().repaint();
        });
        return item;
    }


    public JMenuItem getGenotypeSortItem(final Variant variant) {

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

    public JMenuItem getGenotypeGroupItem() {

        JMenuItem item = new JMenuItem("Group By...");
        item.addActionListener(evt -> (new GroupSamplesMenuAction("", 0, IGV.getInstance())).doGroupBy());
        item.setEnabled(AttributeManager.getInstance().getVisibleAttributes().size() > 0);

        return item;
    }

    public JMenuItem getSampleNameSortItem(final Variant variant) {
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

    public JMenuItem getDepthSortItem(final Variant variant) {
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

    public JMenuItem getQualitySortItem(final Variant variant) {
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

    public void changeVisibilityWindow() {
        TrackMenuUtils.changeFeatureVisibilityWindow(Arrays.asList((Track) track));
    }

    public Collection<JMenuItem> getSortMenuItems(Variant variant) {

        java.util.List<JMenuItem> items = new ArrayList<JMenuItem>();
        items.add(getGenotypeSortItem(variant));
        items.add(getSampleNameSortItem(variant));
        items.add(getDepthSortItem(variant));
        items.add(getQualitySortItem(variant));
        return items;
    }

    public List<JMenuItem> getDisplayModeItems() {

        List<JMenuItem> items = new ArrayList();

        ButtonGroup group = new ButtonGroup();

        Track.DisplayMode displayMode = track.getDisplayMode();

        JRadioButtonMenuItem m1 = new JRadioButtonMenuItem("Collapsed");
        m1.setSelected(displayMode == Track.DisplayMode.COLLAPSED);
        m1.addActionListener(evt -> {
            track.setDisplayMode(Track.DisplayMode.COLLAPSED);
            track.repaint();
        });

        JRadioButtonMenuItem m2 = new JRadioButtonMenuItem("Squished");
        m2.setSelected(displayMode == Track.DisplayMode.SQUISHED);
        m2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.SQUISHED);
                track.repaint();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("Expanded");
        m3.setSelected(displayMode == Track.DisplayMode.EXPANDED);
        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.EXPANDED);
                track.repaint();
            }
        });


        items.add(m1);
        items.add(m2);
        items.add(m3);
        group.add(m1);
        group.add(m2);
        group.add(m3);

        return items;
    }


}
