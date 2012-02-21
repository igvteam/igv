/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.variant;

import org.apache.log4j.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
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

    private static Logger log = Logger.getLogger(VariantMenu.class);
    private VariantTrack track;
    private Collection<String> selectedSamples;
    static boolean depthSortingDirection;
    static boolean genotypeSortingDirection;
    static boolean sampleSortingDirection;
    static boolean qualitySortingDirection;

    public VariantMenu(final VariantTrack variantTrack, Variant variant) {

        this.track = variantTrack;

        if (track.hasAlignmentFiles()) {
            selectedSamples = track.getSelectedSamples();
        }

        this.addPopupMenuListener(new PopupMenuListener() {
            public void popupMenuWillBecomeVisible(PopupMenuEvent popupMenuEvent) {

            }

            public void popupMenuWillBecomeInvisible(PopupMenuEvent popupMenuEvent) {
                close();
            }

            public void popupMenuCanceled(PopupMenuEvent popupMenuEvent) {
                close();
            }

            private void close() {
               // track.clearSelectedVariant();
            }

        });


        //Title
        JLabel popupTitle = new JLabel("<html><b>" + this.track.getName(), JLabel.LEFT);
        Font newFont = getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        add(popupTitle);

        //Change Track Settings
        addSeparator();

        List<Track> selectedTracks = Arrays.asList((Track) variantTrack);
        add(TrackMenuUtils.getTrackRenameItem(selectedTracks));
        add(TrackMenuUtils.getChangeFontSizeItem(selectedTracks));


        //Hides
        addSeparator();
        JLabel colorByItem = new JLabel("<html>&nbsp;&nbsp;<b>Color By", JLabel.LEFT);
        add(colorByItem);
        add(getColorByGenotype());
        add(getColorByAllele());
        if (track.isEnableMethylationRateSupport()) {
            add(getColorByMethylationRate());
        }

        //add(getRenderIDItem());

        //Sorter
        addSeparator();
        for (JMenuItem item : getSortMenuItems(variant)) {
            add(item);
            if (variant == null) {
                item.setEnabled(false);
            }
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
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                int currentValue = track.getSquishedHeight();
                int newValue = TrackMenuUtils.getIntValue("Squished row height", currentValue);
                if (newValue != Integer.MIN_VALUE) {
                    track.setSquishedHeight(newValue);
                    IGV.getInstance().getContentPane().repaint();
                }
            }
        });
        add(item);

        add(getHideFilteredItem());
        add(getFeatureVisibilityItem());

        if (track.hasAlignmentFiles()) {
            addSeparator();
            add(getLoadBamsItem());
        }


        addSeparator();
        add(TrackMenuUtils.getRemoveMenuItem(Arrays.asList(new Track[]{this.track})));


    }

    private JMenuItem getFeatureVisibilityItem() {
        JMenuItem item = new JMenuItem("Set Feature Visibility Window...");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                changeVisibilityWindow();
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }

    private JMenuItem getColorByGenotype() {
        final JMenuItem item = new JCheckBoxMenuItem("Genotype", track.getColorMode() == VariantTrack.ColorMode.GENOTYPE);
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setColorMode(VariantTrack.ColorMode.GENOTYPE);
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }


    private JMenuItem getColorByAllele() {
        final JMenuItem item = new JCheckBoxMenuItem("Allele", track.getColorMode() == VariantTrack.ColorMode.ALLELE);
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setColorMode(VariantTrack.ColorMode.ALLELE);
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }

    private JMenuItem getColorByMethylationRate() {
        final JMenuItem item = new JCheckBoxMenuItem("Methylation Rate", track.getColorMode() == VariantTrack.ColorMode.METHYLATION_RATE);
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setColorMode(VariantTrack.ColorMode.METHYLATION_RATE);
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }


    private JMenuItem getHideFilteredItem() {
        JMenuItem item = new JCheckBoxMenuItem("Suppress Filtered Sites", track.getHideFiltered());
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setHideFiltered(!track.getHideFiltered());
                IGV.getInstance().getContentPane().repaint();
            }
        });
        return item;
    }


    public JMenuItem getGenotypeSortItem(final Variant variant) {

        JMenuItem item = new JMenuItem("Sort By Genotype");
        if (variant != null) {
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    GenotypeComparator compare = new GenotypeComparator(variant);
                    genotypeSortingDirection = !genotypeSortingDirection;
                    track.sortSamples(compare);
                    IGV.getInstance().getContentPane().repaint();
                }
            });
        }

        return item;
    }

    public JMenuItem getSampleNameSortItem(final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Sample Name");
        if (variant != null) {
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
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
                    track.sortSamples(compare);
                    IGV.getInstance().getContentPane().repaint();
                }
            });
        }
        return item;
    }

    public JMenuItem getDepthSortItem(final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Depth");
        if (variant != null) {
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    DepthComparator compare = new DepthComparator(variant);
                    depthSortingDirection = !depthSortingDirection;
                    track.sortSamples(compare);
                    IGV.getInstance().getContentPane().repaint();
                }
            });

        }
        return item;
    }

    public JMenuItem getQualitySortItem(final Variant variant) {
        JMenuItem item = new JMenuItem("Sort By Quality");
        if (variant != null) {
            double quality = variant.getPhredScaledQual();
            if (quality > -1) {
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent evt) {
                        QualityComparator compare = new QualityComparator(variant);
                        qualitySortingDirection = !qualitySortingDirection;
                        track.sortSamples(compare);
                        IGV.getInstance().getContentPane().repaint();
                    }
                });
            } else {
                item.setEnabled(false);
            }
        }

        return item;
    }

    public void changeVisibilityWindow() {
        int value = getIntValue("Visibility Window", track.getVisibilityWindow());
        if (value > 0) {
            track.setVisibilityWindow(value);
        }
    }

    private static int getIntValue(String parameter, int value) {
        while (true) {
            String height = JOptionPane.showInputDialog(
                    IGV.getMainFrame(), parameter + ": ",
                    String.valueOf(value));
            if ((height == null) || height.trim().equals("")) {
                return Integer.MIN_VALUE;   // <= the logical "null" value
            }

            try {
                value = Integer.parseInt(height);
                return value;
            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        parameter + " must be an integer number.");
            }
        }
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
        m1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.COLLAPSED);
                IGV.getInstance().doRefresh();
            }
        });

        JRadioButtonMenuItem m2 = new JRadioButtonMenuItem("Squished");
        m2.setSelected(displayMode == Track.DisplayMode.SQUISHED);
        m2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.SQUISHED);
                IGV.getInstance().doRefresh();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("Expanded");
        m3.setSelected(displayMode == Track.DisplayMode.EXPANDED);
        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.EXPANDED);
                IGV.getInstance().doRefresh();
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

            double depth1 = Double.isNaN(readDepth1) ? -1 : readDepth1;
            double depth2 = Double.isNaN(readDepth2) ? -1 : readDepth2;
            if (depth2 == depth1) {
                return 0;
            } else if (depth2 < depth1) {
                return depthSortingDirection ? -1 : 1;
            } else {
                return depthSortingDirection ? 1 : 1;
            }
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

            if (qual2 == qual1) {
                return 0;
            } else if (qual2 < qual1) {
                return qualitySortingDirection ? -1 : 1;
            } else {
                return qualitySortingDirection ? 1 : 1;
            }
        }
    }


}
