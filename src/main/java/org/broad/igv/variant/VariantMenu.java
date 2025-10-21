/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.variant;

import org.broad.igv.logging.*;
import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.GroupTracksMenuAction;
import org.broad.igv.ui.panel.IGVPopupMenu;

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
    private Collection<String> selectedSamples;
    static boolean depthSortingDirection;
    static boolean genotypeSortingDirection;
    static boolean sampleSortingDirection;
    static boolean qualitySortingDirection;


    public VariantMenu(final VariantTrack variantTrack, final Variant variant, TrackClickEvent e) {
        super();
        this.track = variantTrack;

        if (track.hasAlignmentFiles()) {
            selectedSamples = track.getSelectedSamples();
        }

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
            addSeparator();
            add(getGenotypeGroupItem());
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

        if (track.hasAlignmentFiles()) {
            addSeparator();
            add(getLoadBamsItem());
        }
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
                GenotypeComparator compare = new GenotypeComparator(variant);
                genotypeSortingDirection = !genotypeSortingDirection;
                track.sortSamplesByAttribute(compare);
                IGV.getInstance().getContentPane().repaint();
            });
        }

        return item;
    }

    public JMenuItem getGenotypeGroupItem() {

        JMenuItem item = new JMenuItem("Group By...");
        item.addActionListener(evt -> (new GroupTracksMenuAction("", 0, IGV.getInstance())).doGroupBy());
        item.setEnabled(AttributeManager.getInstance().getVisibleAttributes().size() > 0);

        return item;
    }

    public JMenuItem getSampleNameSortItem(final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Sample Name");
        if (variant != null) {
            item.addActionListener(evt -> {
                Comparator<String> compare = new Comparator<String>() {
                    public int compare(String o, String o1) {
                        if (sampleSortingDirection) {
                            return o.compareTo(o1);
                        } else {
                            return o1.compareTo(o);
                        }
                    }
                };
                sampleSortingDirection = !sampleSortingDirection;
                track.sortSamplesByAttribute(compare);
                IGV.getInstance().getContentPane().repaint();
            });
        }
        return item;
    }

    public JMenuItem getDepthSortItem(final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Depth");
        if (variant != null) {
            item.addActionListener(evt -> {
                DepthComparator compare = new DepthComparator(variant);
                depthSortingDirection = !depthSortingDirection;
                track.sortSamplesByAttribute(compare);
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
                    QualityComparator compare = new QualityComparator(variant);
                    qualitySortingDirection = !qualitySortingDirection;
                    track.sortSamplesByAttribute(compare);
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

    /**
     * Load bam files associated with the selected samples (experimental).
     *
     * @return
     */
    private JMenuItem getLoadBamsItem() {
        final JMenuItem item = new JMenuItem("Load alignments");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.loadSelectedBams();
            }
        });
        item.setEnabled(selectedSamples != null && selectedSamples.size() > 0);
        return item;
    }


    static class GenotypeComparator implements Comparator<String> {

        Variant variant;

        GenotypeComparator(Variant variant) {
            this.variant = variant;
        }

        public int compare(String e1, String e2) {

            int genotype1 = classifyGenotype(variant.getGenotype(e1));
            int genotype2 = classifyGenotype(variant.getGenotype(e2));

            if (genotype2 == genotype1) {
                return 0;
            } else if (genotype2 > genotype1) {
                return genotypeSortingDirection ? 1 : -1;
            } else {
                return genotypeSortingDirection ? -1 : 1;
            }
        }


        private int classifyGenotype(Genotype genotype) {

            if (genotype.isNoCall()) {
                return genotypeSortingDirection ? 1 : 10;
            } else if (genotype.isHomVar()) {
                return 4;
            } else if (genotype.isHet()) {
                return 3;
            } else if (genotype.isHomRef()) {
                return genotypeSortingDirection ? 2 : 9;
            }
            return -1; //Unknown
        }
    }


    static class DepthComparator implements Comparator<String> {

        Variant variant;

        DepthComparator(Variant variant) {
            this.variant = variant;
        }

        public int compare(String s1, String s2) {


            double readDepth1 = variant.getGenotype(s1).getAttributeAsDouble("DP");
            double readDepth2 = variant.getGenotype(s2).getAttributeAsDouble("DP");

            int sign = depthSortingDirection ? -1 : 1;
            return sign * Double.compare(readDepth1, readDepth2);

        }
    }

    static class QualityComparator implements Comparator<String> {

        Variant variant;

        QualityComparator(Variant variant) {
            this.variant = variant;
        }

        public int compare(String s1, String s2) {

            double qual1 = variant.getGenotype(s1).getPhredScaledQual();
            double qual2 = variant.getGenotype(s2).getPhredScaledQual();

            int sign = qualitySortingDirection ? -1 : 1;
            return sign * Double.compare(qual1, qual2);

        }
    }


}
