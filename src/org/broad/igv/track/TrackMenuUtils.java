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

package org.broad.igv.track;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.CombinedDataSource;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.basepair.BasePairTrack;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.renderer.*;
import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.sam.SAMWriter;
import org.broad.igv.ui.*;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.collections.CollUtils;
import org.broad.igv.util.stats.KMPlotFrame;
import htsjdk.tribble.Feature;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class TrackMenuUtils {

    static Logger log = Logger.getLogger(TrackMenuUtils.class);
    final static String LEADING_HEADING_SPACER = "  ";
    private static final WindowFunction[] ORDERED_WINDOW_FUNCTIONS = new WindowFunction[]{
            WindowFunction.min,
            WindowFunction.percentile2,
            WindowFunction.percentile10,
            WindowFunction.median,
            WindowFunction.mean,
            WindowFunction.percentile90,
            WindowFunction.percentile98,
            WindowFunction.max,
            WindowFunction.none
    };

    private static List<TrackMenuItemBuilder> trackMenuItems = new ArrayList<TrackMenuItemBuilder>();

    /**
     * Called by plugins to add a listener, which is then called when TrackMenus are created
     * to generate menu entries.
     *
     * @param builder
     * @api
     */
    public static void addTrackMenuItemBuilder(TrackMenuItemBuilder builder) {
        trackMenuItems.add(builder);
    }


    /**
     * Return a popup menu with items applicable to the collection of tracks.
     *
     * @param tracks
     * @return
     */
    public static IGVPopupMenu getPopupMenu(final Collection<Track> tracks, String title, TrackClickEvent te) {

        if (log.isDebugEnabled()) {
            log.debug("enter getPopupMenu");
        }

        IGVPopupMenu menu = new IGVPopupMenu();

        JLabel popupTitle = new JLabel(LEADING_HEADING_SPACER + title, JLabel.CENTER);
        popupTitle.setFont(UIConstants.boldFont);
        if (popupTitle != null) {
            menu.add(popupTitle);
            menu.addSeparator();
        }

        addStandardItems(menu, tracks, te);

        return menu;

    }

    /**
     * Add menu items which have been added through the api, not known until runtime
     *
     * @param menu
     * @param tracks
     * @param te
     */
    public static void addPluginItems(JPopupMenu menu, Collection<Track> tracks, TrackClickEvent te) {
        List<JMenuItem> items = new ArrayList<JMenuItem>(0);
        for (TrackMenuItemBuilder builder : trackMenuItems) {
            JMenuItem item = builder.build(tracks, te);
            if (item != null) {
                items.add(item);
            }
        }

        if (items.size() > 0) {
            menu.addSeparator();
            for (JMenuItem item : items) {
                menu.add(item);
            }
        }
    }

    public static void addStandardItems(JPopupMenu menu, Collection<Track> tracks, TrackClickEvent te) {

        boolean hasDataTracks = false;
        boolean hasFeatureTracks = false;
        boolean hasOtherTracks = false;
        boolean hasCoverageTracks = false;
        for (Track track : tracks) {

            if (track instanceof DataTrack) {
                hasDataTracks = true;
            } else if(track instanceof CoverageTrack) {
                hasDataTracks = true;
                hasCoverageTracks = true;
            } else if (track instanceof FeatureTrack) {
                hasFeatureTracks = true;
            } else {
                hasOtherTracks = true;
            }
            if (hasDataTracks && hasFeatureTracks && hasOtherTracks) {
                break;
            }
        }

        boolean hasBasePairTracks = false;
        for (Track track : tracks) {
            if (track instanceof BasePairTrack) {
                hasBasePairTracks = true;
                break;
            }
        }
        if (hasBasePairTracks) {
            addBasePairItems(menu, tracks);
        }

        boolean featureTracksOnly = hasFeatureTracks && !hasDataTracks && !hasOtherTracks;
        boolean dataTracksOnly = !hasFeatureTracks && hasDataTracks && !hasOtherTracks;

        addSharedItems(menu, tracks, hasFeatureTracks, hasCoverageTracks);
        menu.addSeparator();
        if (dataTracksOnly) {
            addDataItems(menu, tracks, hasCoverageTracks);
        } else if (featureTracksOnly) {
            addFeatureItems(menu, tracks, te);
        }


    }

    public static void addZoomItems(JPopupMenu menu, final ReferenceFrame frame) {

        if (FrameManager.isGeneListMode()) {
            JMenuItem item = new JMenuItem("Reset panel to '" + frame.getName() + "'");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    frame.reset();
                    // TODO -- paint only panels for this frame
                }
            });
            menu.add(item);
        }


        JMenuItem zoomOutItem = new JMenuItem("Zoom out");
        zoomOutItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                frame.doZoomIncrement(-1);
            }
        });
        menu.add(zoomOutItem);

        JMenuItem zoomInItem = new JMenuItem("Zoom in");
        zoomInItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                frame.doZoomIncrement(1);
            }
        });
        menu.add(zoomInItem);


    }


    /**
     * Return popup menu with items applicable to data tracks
     *
     * @return
     */
    public static void addDataItems(JPopupMenu menu, final Collection<Track> tracks, boolean hasCoverageTracks) {

        if (log.isTraceEnabled()) {
            log.trace("enter getDataPopupMenu");
        }

        if(!hasCoverageTracks) {
            final String[] labels = {"Heatmap", "Bar Chart", "Points", "Line Plot"};
            final Class[] renderers = {HeatmapRenderer.class, BarChartRenderer.class,
                    PointsRenderer.class, LineplotRenderer.class
            };

            JLabel rendererHeading = new JLabel(LEADING_HEADING_SPACER + "Type of Graph", JLabel.LEFT);
            rendererHeading.setFont(UIConstants.boldFont);

            menu.add(rendererHeading);

            // Get existing selections
            Set<Class> currentRenderers = new HashSet<Class>();
            for (Track track : tracks) {
                if (track.getRenderer() != null) {
                    currentRenderers.add(track.getRenderer().getClass());
                }
            }

            // Create renderer menu items
            for (int i = 0; i < labels.length; i++) {
                JCheckBoxMenuItem item = new JCheckBoxMenuItem(labels[i]);
                final Class rendererClass = renderers[i];
                if (currentRenderers.contains(rendererClass)) {
                    item.setSelected(true);
                }
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent evt) {
                        changeRenderer(tracks, rendererClass);
                    }
                });
                menu.add(item);
            }
            menu.addSeparator();


            // Get union of all valid window functions for selected tracks
            Set<WindowFunction> avaibleWindowFunctions = new HashSet();
            for (Track track : tracks) {
                avaibleWindowFunctions.addAll(track.getAvailableWindowFunctions());
            }
            avaibleWindowFunctions.add(WindowFunction.none);


            // dataPopupMenu.addSeparator();
            // Collection all window functions for selected tracks
            Set<WindowFunction> currentWindowFunctions = new HashSet<WindowFunction>();
            for (Track track : tracks) {
                if (track.getWindowFunction() != null) {
                    currentWindowFunctions.add(track.getWindowFunction());
                }
            }

            if (avaibleWindowFunctions.size() > 1 || currentWindowFunctions.size() > 1) {
                JLabel statisticsHeading = new JLabel(LEADING_HEADING_SPACER + "Windowing Function", JLabel.LEFT);
                statisticsHeading.setFont(UIConstants.boldFont);

                menu.add(statisticsHeading);

                for (final WindowFunction wf : ORDERED_WINDOW_FUNCTIONS) {
                    JCheckBoxMenuItem item = new JCheckBoxMenuItem(wf.getValue());
                    if (avaibleWindowFunctions.contains(wf) || currentWindowFunctions.contains(wf)) {
                        if (currentWindowFunctions.contains(wf)) {
                            item.setSelected(true);
                        }
                        item.addActionListener(new ActionListener() {

                            public void actionPerformed(ActionEvent evt) {
                                changeStatType(wf.toString(), tracks);
                            }
                        });
                        menu.add(item);
                    }
                }
                menu.addSeparator();
            }
        }


        menu.add(getDataRangeItem(tracks));

        if(!hasCoverageTracks) menu.add(getHeatmapScaleItem(tracks));

        if (tracks.size() > 0) {
            menu.add(getLogScaleItem(tracks));
        }

        menu.add(getAutoscaleItem(tracks));

        if (tracks.size() > 1 ) {
            menu.add(getGroupAutoscaleItem(tracks));
        }

        menu.add(getShowDataRangeItem(tracks));

        //Optionally add overlay track options
        Track firstTrack = tracks.iterator().next();
        boolean merged = (tracks.size() == 1 && firstTrack instanceof MergedTracks);

        if(tracks.size() > 1 || merged) {
            menu.addSeparator();
            final List<DataTrack> dataTrackList = Lists.newArrayList(Iterables.filter(tracks, DataTrack.class));
            final JMenuItem overlayGroups = new JMenuItem("Overlay Tracks");
            overlayGroups.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    MergedTracks mergedTracks = new MergedTracks(UUID.randomUUID().toString(), "Overlay", dataTrackList);
                    Track firstTrack = tracks.iterator().next();
                    TrackPanel panel = TrackPanel.getParentPanel(firstTrack);
                    panel.addTrack(mergedTracks);
                    panel.moveSelectedTracksTo(Arrays.asList(mergedTracks), firstTrack, false);
                    panel.removeTracks(tracks);
                }
            });

            int numDataTracks = dataTrackList.size();
            overlayGroups.setEnabled(numDataTracks >= 2 && numDataTracks == tracks.size());
            menu.add(overlayGroups);

            // Enable "separateTracks" menu if selection is a single track, and that track is merged.

            JMenuItem unmergeItem = new JMenuItem("Separate Tracks");
            menu.add(unmergeItem);

            if (merged) {

                unmergeItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        Track firstTrack = tracks.iterator().next();
                        TrackPanel panel = TrackPanel.getParentPanel(firstTrack);
                        final MergedTracks mergedTracks = (MergedTracks) firstTrack;
                        mergedTracks.setTrackAlphas(255);
                        panel.addTracks(mergedTracks.getMemberTracks());
                        panel.moveSelectedTracksTo(mergedTracks.getMemberTracks(), mergedTracks, true);
                        IGV.getInstance().removeTracks(Arrays.asList(mergedTracks));
                    }
                });
            } else {
                unmergeItem.setEnabled(false);
            }
        }

        //menu.addSeparator();
        //menu.add(getChangeKMPlotItem(tracks));

        if (Globals.isDevelopment() && FrameManager.isGeneListMode() && tracks.size() == 1) {
            menu.addSeparator();
            menu.add(getShowSortFramesItem(tracks.iterator().next()));
        }


    }

    private static List<JMenuItem> getCombinedDataSourceItems(final Collection<Track> tracks) {

        Iterable<DataTrack> dataTracksIter = Iterables.filter(tracks, DataTrack.class);
        final List<DataTrack> dataTracks = Lists.newArrayList(dataTracksIter);
        JMenuItem addItem = new JMenuItem("Sum Tracks");
        JMenuItem subItem = new JMenuItem("Subtract Tracks");
        boolean enableComb = dataTracks.size() == 2;

        addItem.setEnabled(enableComb);
        addItem.setEnabled(enableComb);

        addItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                addCombinedDataTrack(dataTracks, CombinedDataSource.Operation.ADD);
            }
        });

        subItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                addCombinedDataTrack(dataTracks, CombinedDataSource.Operation.SUBTRACT);
            }
        });

        return Arrays.asList(addItem, subItem);
    }

    private static void addCombinedDataTrack(List<DataTrack> dataTracks, CombinedDataSource.Operation op) {
        String text = "";
        switch (op) {
            case ADD:
                text = "Sum";
                break;
            case SUBTRACT:
                text = "Difference";
                break;
        }
        DataTrack track0 = dataTracks.get(0);
        DataTrack track1 = dataTracks.get(1);
        CombinedDataSource source = new CombinedDataSource(track0, track1, op);

        DataSourceTrack newTrack = new DataSourceTrack(null, track0.getId() + track1.getId() + text, text, source);
        changeRenderer(Arrays.<Track>asList(newTrack), track0.getRenderer().getClass());
        newTrack.setDataRange(track0.getDataRange());
        newTrack.setColorScale(track0.getColorScale());
        IGV.getInstance().addTracks(Arrays.<Track>asList(newTrack), PanelName.DATA_PANEL);
    }

    /**
     * Return popup menu with items applicable to feature tracks
     *
     * @return
     */
    private static void addFeatureItems(JPopupMenu featurePopupMenu, final Collection<Track> tracks, TrackClickEvent te) {


        addDisplayModeItems(tracks, featurePopupMenu);


        if (tracks.size() == 1) {
            Track t = tracks.iterator().next();
            Feature f = t.getFeatureAtMousePosition(te);
            if (f != null) {
                featurePopupMenu.addSeparator();

                // If we are over an exon, copy its sequence instead of the entire feature.
                if (f instanceof IGVFeature) {
                    double position = te.getChromosomePosition();
                    Collection<Exon> exons = ((IGVFeature) f).getExons();
                    if (exons != null) {
                        for (Exon exon : exons) {
                            if (position > exon.getStart() && position < exon.getEnd()) {
                                f = exon;
                                break;
                            }
                        }
                    }
                }


                featurePopupMenu.add(getCopyDetailsItem(f, te));
                featurePopupMenu.add(getCopySequenceItem(f));
                featurePopupMenu.add(getBlatItem(f));


            }
            if (Globals.isDevelopment()) {
                featurePopupMenu.addSeparator();
                featurePopupMenu.add(getFeatureToGeneListItem(t));
            }
            if (Globals.isDevelopment() && FrameManager.isGeneListMode() && tracks.size() == 1) {
                featurePopupMenu.addSeparator();
                featurePopupMenu.add(getShowSortFramesItem(tracks.iterator().next()));
            }

        }

        featurePopupMenu.addSeparator();
        featurePopupMenu.add(getChangeFeatureWindow(tracks));

    }

    /**
     * Return popup menu with items applicable to arc tracks
     *
     * @return
     */
    // stevenbusan
    public static void addBasePairItems(JPopupMenu menu, final Collection<Track> tracks) {

        JLabel arcDirectionHeading = new JLabel(LEADING_HEADING_SPACER + "Arc direction", JLabel.LEFT);
        arcDirectionHeading.setFont(UIConstants.boldFont);

        menu.add(arcDirectionHeading);

        final String[] arcDirectionLabels = {"Up", "Down"};

        for (int i = 0; i < arcDirectionLabels.length; i++) {
            JCheckBoxMenuItem item = new JCheckBoxMenuItem(arcDirectionLabels[i]);
            final int n = (i==0) ? 1 : -1 ;
            for (Track track : tracks) {
                if (track instanceof BasePairTrack) {
                    if (((BasePairTrack) track).getDirection() == n){
                        item.setSelected(true);
                    }
                }
            }
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    for (Track track : tracks) {
                        if (track instanceof BasePairTrack) {
                            ((BasePairTrack) track).setDirection(n);
                        }
                    }
                    IGV.getInstance().repaint();
                }
            });
            menu.add(item);
        }

        menu.addSeparator();
    }

    private static JMenuItem getFeatureToGeneListItem(final Track t) {
        JMenuItem mi = new JMenuItem("Use as loci list");

        mi.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Current chromosome only for now


            }
        });

        return mi;
    }

    /**
     * Return a menu item which will export visible features
     * If {@code tracks} is not a single {@code FeatureTrack}, {@code null}
     * is returned (there should be no menu entry)
     *
     * @param tracks
     * @return
     */
    public static JMenuItem getExportFeatures(final Collection<Track> tracks, final ReferenceFrame frame) {
        Track ft = tracks.iterator().next();
        if (tracks.size() != 1) {
            return null;
        }
        JMenuItem exportData = null;

        if (ft instanceof FeatureTrack) {
            exportData = new JMenuItem("Export Features...");
            exportData.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    File outFile = FileDialogUtils.chooseFile("Save Visible Data",
                            PreferenceManager.getInstance().getLastTrackDirectory(),
                            new File("visibleData.bed"),
                            FileDialogUtils.SAVE);

                    exportVisibleFeatures(outFile.getAbsolutePath(), tracks, frame);
                }
            });
        } else if (ft instanceof AlignmentTrack) {
            exportData = new JMenuItem("Export Alignments...");
            exportData.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    File outFile = FileDialogUtils.chooseFile("Save Visible Data",
                            PreferenceManager.getInstance().getLastTrackDirectory(),
                            new File("visibleData.sam"),
                            FileDialogUtils.SAVE);

                    int countExp = exportVisibleAlignments(outFile.getAbsolutePath(), tracks, frame);
                    String msg = String.format("%d reads written", countExp);
                    MessageUtils.setStatusBarMessage(msg);
                }
            });
        }

        return exportData;
    }

    static int exportVisibleAlignments(String outPath, Collection<Track> tracks, ReferenceFrame frame) {
        AlignmentTrack alignmentTrack = null;
        for (Track track : tracks) {
            if (track instanceof AlignmentTrack) {
                alignmentTrack = (AlignmentTrack) track;
                break;
            }
        }

        if (alignmentTrack == null) return -1;


        File outFile = new File(outPath);
        try {
            AlignmentDataManager dataManager = alignmentTrack.getDataManager();
            ResourceLocator inlocator = dataManager.getLocator();
            Range range = frame.getCurrentRange();

            //Read directly from file
            //return SAMWriter.writeAlignmentFilePicard(inlocator, outPath, range.getChr(), range.getStart(), range.getEnd());

            //Export those in memory, overlapping current view
            return SAMWriter.writeAlignmentFilePicard(dataManager, outFile, range.getChr(), range.getStart(), range.getEnd());
        } catch (IOException e) {
            log.error(e.getMessage(), e);
            throw new RuntimeException(e);
        }

    }

    /**
     * Write features in {@code track} found in {@code range} to {@code outPath},
     * BED format
     * TODO Move somewhere else? run on separate thread?  Probably shouldn't be here
     *
     * @param outPath
     * @param tracks
     * @param frame
     */
    static void exportVisibleFeatures(String outPath, Collection<Track> tracks, ReferenceFrame frame) {
        PrintWriter writer;
        try {
            writer = new PrintWriter(outPath);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }

        for (Track track : tracks) {
            if (track instanceof FeatureTrack) {
                FeatureTrack fTrack = (FeatureTrack) track;

                String trackLine = fTrack.getExportTrackLine();
                if(trackLine != null) {
                    writer.println(trackLine);
                }

                //Can't trust FeatureTrack.getFeatures to limit itself, so we filter
                List<Feature> features = fTrack.getVisibleFeatures(frame);
                Range range = frame.getCurrentRange();
                Predicate<Feature> pred = FeatureUtils.getOverlapPredicate(range.getChr(), range.getStart(), range.getEnd());
                features = CollUtils.filter(features, pred);
                IGVBEDCodec codec = new IGVBEDCodec();
                for (Feature feat : features) {
                    String featString = codec.encode(feat);
                    writer.println(featString);
                }
            }
        }
        writer.flush();
        writer.close();
    }


    /**
     * Popup menu with items applicable to both feature and data tracks
     *
     * @return
     */
    public static void addSharedItems(JPopupMenu menu, final Collection<Track> tracks, boolean hasFeatureTracks, boolean hasCoverageTracks) {

        //JLabel trackSettingsHeading = new JLabel(LEADING_HEADING_SPACER + "Track Settings", JLabel.LEFT);
        //trackSettingsHeading.setFont(boldFont);
        //menu.add(trackSettingsHeading);

        menu.add(getTrackRenameItem(tracks));

        String colorLabel = (hasFeatureTracks || hasCoverageTracks)
                ? "Change Track Color..." : "Change Track Color (Positive Values)...";
        JMenuItem item = new JMenuItem(colorLabel);
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                changeTrackColor(tracks);
            }
        });
        menu.add(item);

        if (!(hasFeatureTracks || hasCoverageTracks)) {

            // Change track color by attribute
            item = new JMenuItem("Change Track Color (Negative Values)...");
            item.setToolTipText(
                    "Change the alternate track color.  This color is used when graphing negative values");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    changeAltTrackColor(tracks);
                }
            });
            menu.add(item);
        }

        menu.add(getChangeTrackHeightItem(tracks));
        menu.add(getChangeFontSizeItem(tracks));
    }


    private static void changeStatType(String statType, Collection<Track> selectedTracks) {
        for (Track track : selectedTracks) {
            track.setWindowFunction(WindowFunction.valueOf(statType));
        }
        refresh();
    }


    public static JMenuItem getTrackRenameItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Rename Track...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        renameTrack(selectedTracks);
                    }
                });
            }
        });
        if (selectedTracks.size() > 1) {
            item.setEnabled(false);
        }
        return item;
    }

    private static JMenuItem getHeatmapScaleItem(final Collection<Track> selectedTracks) {

        JMenuItem item = new JMenuItem("Set Heatmap Scale...");

        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                if (selectedTracks.size() > 0) {

                    ContinuousColorScale colorScale = selectedTracks.iterator().next().getColorScale();
                    HeatmapScaleDialog dlg = new HeatmapScaleDialog(IGV.getMainFrame(), colorScale);

                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        colorScale = dlg.getColorScale();

                        // dlg.isFlipAxis());
                        for (Track track : selectedTracks) {
                            track.setColorScale(colorScale);
                        }
                        IGV.getInstance().repaint();
                    }

                }

            }
        });
        return item;
    }

    public static JMenuItem getDataRangeItem(final Collection<Track> selectedTracks) {
        JMenuItem item = new JMenuItem("Set Data Range...");

        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                if (selectedTracks.size() > 0) {

                    // Create a datarange that spans the extent of prev tracks range
                    DataRange prevAxisDefinition = DataRange.getFromTracks(selectedTracks);
                    DataRangeDialog dlg = new DataRangeDialog(IGV.getMainFrame(), prevAxisDefinition);
                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        float min = Math.min(dlg.getMax(), dlg.getMin());
                        float max = Math.max(dlg.getMin(), dlg.getMax());
                        float mid = dlg.getBase();
                        mid = Math.max(min, Math.min(mid, max));

                        DataRange axisDefinition = new DataRange(dlg.getMin(), mid, dlg.getMax(),
                                prevAxisDefinition.isDrawBaseline(), dlg.isLog());

                        for (Track track : selectedTracks) {
                            track.setDataRange(axisDefinition);
                            track.setAutoScale(false);
                            track.removeAttribute(AttributeManager.GROUP_AUTOSCALE);
                        }
                        IGV.getInstance().repaint();
                    }

                }
            }
        });

        return item;
    }

    private static JMenuItem getDrawBorderItem() {
        // Change track height by attribute


        final JCheckBoxMenuItem drawBorderItem = new JCheckBoxMenuItem("Draw borders");
        drawBorderItem.setSelected(FeatureTrack.isDrawBorder());
        drawBorderItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                FeatureTrack.setDrawBorder(drawBorderItem.isSelected());
                IGV.getInstance().repaintDataPanels();
            }
        });

        return drawBorderItem;
    }


    public static JMenuItem getLogScaleItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute


        final JCheckBoxMenuItem logScaleItem = new JCheckBoxMenuItem("Log scale");
        final boolean logScale = selectedTracks.iterator().next().getDataRange().isLog();
        logScaleItem.setSelected(logScale);
        logScaleItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                DataRange.Type scaleType = logScaleItem.isSelected() ?
                        DataRange.Type.LOG :
                        DataRange.Type.LINEAR;
                for (Track t : selectedTracks) {
                    t.getDataRange().setType(scaleType);
                }
                IGV.getInstance().repaintDataPanels();
            }
        });

        return logScaleItem;
    }

    private static JMenuItem getAutoscaleItem(final Collection<Track> selectedTracks) {

        final JCheckBoxMenuItem autoscaleItem = new JCheckBoxMenuItem("Autoscale");
        if (selectedTracks.size() == 0) {
            autoscaleItem.setEnabled(false);

        } else {
            boolean autoScale = checkAutoscale(selectedTracks);

            autoscaleItem.setSelected(autoScale);
            autoscaleItem.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent evt) {

                    boolean autoScale = autoscaleItem.isSelected();
                    for (Track t : selectedTracks) {
                        t.setAutoScale(autoScale);
                        if (autoScale) {
                            t.removeAttribute(AttributeManager.GROUP_AUTOSCALE);
                        }
                    }

                    DataPanelContainer.resetStateHash();

                    IGV.getInstance().repaint();
                }
            });
        }
        return autoscaleItem;
    }

    private static JMenuItem getGroupAutoscaleItem(final Collection<Track> selectedTracks) {

        final JMenuItem autoscaleItem = new JMenuItem("Group Autoscale");

        autoscaleItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                int nextAutoscaleGroup = IGV.getInstance().getSession().getNextAutoscaleGroup();
                for (Track t : selectedTracks) {
                    t.setAttributeValue(AttributeManager.GROUP_AUTOSCALE, String.valueOf(nextAutoscaleGroup));
                    t.setAutoScale(false);
                }
                IGV.getInstance().getSession().incrementNextAutoscaleGroup();

                DataPanelContainer.resetStateHash();

                PreferenceManager.getInstance().setShowAttributeView(true);
                IGV.getInstance().getMainPanel().invalidate();
                IGV.getInstance().doRefresh();


            }
        });

        return autoscaleItem;
    }

    private static boolean checkAutoscale(Collection<Track> selectedTracks) {
        boolean autoScale = false;
        for (Track t : selectedTracks) {
            if (t.getAutoScale()) {
                autoScale = true;
                break;
            }
        }
        return autoScale;
    }

    public static JMenuItem getShowDataRangeItem(final Collection<Track> selectedTracks) {

        final JCheckBoxMenuItem item = new JCheckBoxMenuItem("Show Data Range");
        if (selectedTracks.size() == 0) {
            item.setEnabled(false);

        } else {
            boolean showDataRange = true;
            for (Track t : selectedTracks) {
                if (!t.isShowDataRange()) {
                    showDataRange = false;
                    break;
                }
            }

            item.setSelected(showDataRange);
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent evt) {

                    boolean showDataRange = item.isSelected();
                    for (Track t : selectedTracks) {
                        if (t instanceof DataTrack) {
                            ((DataTrack) t).setShowDataRange(showDataRange);
                        }
                    }
                    IGV.getInstance().repaintDataPanels();
                }
            });
        }
        return item;
    }


    public static void addDisplayModeItems(final Collection<Track> tracks, JPopupMenu menu) {

        // Find "most representative" state from track collection
        Map<Track.DisplayMode, Integer> counts = new HashMap<Track.DisplayMode, Integer>(Track.DisplayMode.values().length);
        Track.DisplayMode currentMode = null;

        for (Track t : tracks) {
            Track.DisplayMode mode = t.getDisplayMode();
            if (counts.containsKey(mode)) {
                counts.put(mode, counts.get(mode) + 1);
            } else {
                counts.put(mode, 1);
            }
        }

        int maxCount = -1;
        for (Map.Entry<Track.DisplayMode, Integer> count : counts.entrySet()) {
            if (count.getValue() > maxCount) {
                currentMode = count.getKey();
                maxCount = count.getValue();
            }
        }

        ButtonGroup group = new ButtonGroup();
        Map<String, Track.DisplayMode> modes = new LinkedHashMap<String, Track.DisplayMode>(4);
        modes.put("Collapsed", Track.DisplayMode.COLLAPSED);
        modes.put("Expanded", Track.DisplayMode.EXPANDED);
        modes.put("Squished", Track.DisplayMode.SQUISHED);

        for (final Map.Entry<String, Track.DisplayMode> entry : modes.entrySet()) {
            JRadioButtonMenuItem mm = new JRadioButtonMenuItem(entry.getKey());
            mm.setSelected(currentMode == entry.getValue());
            mm.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    setTrackDisplayMode(tracks, entry.getValue());
                    refresh();
                }
            });
            group.add(mm);
            menu.add(mm);
        }

    }


    private static void setTrackDisplayMode(Collection<Track> tracks, Track.DisplayMode mode) {

        for (Track t : tracks) {
            t.setDisplayMode(mode);
        }
    }


    public static JMenuItem getRemoveMenuItem(final Collection<Track> selectedTracks) {

        boolean multiple = selectedTracks.size() > 1;

        JMenuItem item = new JMenuItem("Remove Track" + (multiple ? "s" : ""));
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                removeTracksAction(selectedTracks);
            }
        });
        return item;
    }

    /**
     * Display a dialog to the user asking to confirm if they want to remove the
     * selected tracks
     *
     * @param selectedTracks
     */
    public static void removeTracksAction(final Collection<Track> selectedTracks) {
        if (selectedTracks.isEmpty()) {
            return;
        }

        StringBuffer buffer = new StringBuffer();
        for (Track track : selectedTracks) {
            buffer.append("\n\t");
            buffer.append(track.getName());
        }
        String deleteItems = buffer.toString();

        JTextArea textArea = new JTextArea();
        textArea.setEditable(false);
        JScrollPane scrollPane = new JScrollPane(textArea);
        textArea.setText(deleteItems);

        JOptionPane optionPane = new JOptionPane(scrollPane,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.YES_NO_OPTION);
        optionPane.setPreferredSize(new Dimension(550, 500));
        JDialog dialog = optionPane.createDialog(IGV.getMainFrame(), "Remove The Following Tracks");
        dialog.setVisible(true);

        Object choice = optionPane.getValue();
        if ((choice == null) || (JOptionPane.YES_OPTION != ((Integer) choice).intValue())) {
            return;
        }

        DataPanelContainer.resetStateHash();
        IGV.getInstance().removeTracks(selectedTracks);
        IGV.getInstance().doRefresh();
    }


    public static void changeRenderer(final Collection<Track> selectedTracks, Class rendererClass) {
        for (Track track : selectedTracks) {

            // TODO -- a temporary hack to facilitate RNAi development
            if (track.getTrackType() == TrackType.RNAI) {
                if (rendererClass == BarChartRenderer.class) {
                    rendererClass = RNAiBarChartRenderer.class;
                }

            }
            track.setRendererClass(rendererClass);
        }
        refresh();
    }

    public static void renameTrack(final Collection<Track> selectedTracks) {

        if (selectedTracks.isEmpty()) {
            return;
        }
        Track t = selectedTracks.iterator().next();
        String newName = JOptionPane.showInputDialog(IGV.getMainFrame(), "Enter new name: ", t.getName());

        if (newName == null || newName.trim() == "") {
            return;
        }

        t.setName(newName);
        refresh();
    }

    public static void changeTrackHeight(final Collection<Track> selectedTracks) {
        if (selectedTracks.isEmpty()) {
            return;
        }

        final String parameter = "Track height";
        Integer value = getIntegerInput(parameter, getRepresentativeTrackHeight(selectedTracks));
        if (value == null) {
            return;
        }

        value = Math.max(0, value);
        for (Track track : selectedTracks) {
            track.setHeight(value, true);
        }
        refresh();
    }

    public static void changeFeatureVisibilityWindow(final Collection<Track> selectedTracks) {

        Collection<Track> featureTracks = new ArrayList(selectedTracks.size());
        for (Track t : selectedTracks) {
            if (t instanceof FeatureTrack) {
                featureTracks.add(t);
            }
        }

        if (featureTracks.isEmpty()) {
            return;
        }


        int origValue = featureTracks.iterator().next().getVisibilityWindow();
        double origValueKB = (origValue / 1000.0);
        Double value = getDoubleInput("Enter visibility window in kilo-bases.  To load all data enter zero.", origValueKB);
        if (value == null) {
            return;
        }

        for (Track track : featureTracks) {
            track.setVisibilityWindow((int) (value * 1000));
        }

        refresh();
    }

    public static void changeFontSize(final Collection<Track> selectedTracks) {


        if (selectedTracks.isEmpty()) {
            return;
        }

        final String parameter = "Font size";
        int defaultValue = selectedTracks.iterator().next().getFontSize();
        Integer value = getIntegerInput(parameter, defaultValue);
        if (value == null) {
            return;
        }

        for (Track track : selectedTracks) {
            track.setFontSize(value);
        }

        refresh();
    }


    public static Integer getIntegerInput(String parameter, int value) {

        while (true) {

            String strValue = JOptionPane.showInputDialog(
                    IGV.getMainFrame(), parameter + ": ",
                    String.valueOf(value));

            //strValue will be null if dialog cancelled
            if ((strValue == null) || strValue.trim().equals("")) {
                return null;
            }

            try {
                value = Integer.parseInt(strValue);
                return value;
            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        parameter + " must be an integer number.");
            }
        }
    }

    public static Double getDoubleInput(String parameter, double value) {

        while (true) {

            String strValue = JOptionPane.showInputDialog(
                    IGV.getMainFrame(), parameter + ": ",
                    String.valueOf(value));

            //strValue will be null if dialog cancelled
            if ((strValue == null) || strValue.trim().equals("")) {
                return null;
            }

            try {
                value = Double.parseDouble(strValue);
                return value;
            } catch (NumberFormatException numberFormatException) {
                MessageUtils.showMessage(parameter + " must be a number.");
            }
        }
    }

    public static void changeTrackColor(final Collection<Track> selectedTracks) {

        if (selectedTracks.isEmpty()) {
            return;
        }

        Color currentSelection = selectedTracks.iterator().next().getColor();

        Color color = UIUtilities.showColorChooserDialog(
                "Select Track Color (Positive Values)",
                currentSelection);

        if (color == null) {
            return;
        }

        for (Track track : selectedTracks) {
            //We preserve the alpha value. This is motivated by MergedTracks
            track.setColor(ColorUtilities.modifyAlpha(color, currentSelection.getAlpha()));
        }
        refresh();

    }

    public static void changeAltTrackColor(final Collection<Track> selectedTracks) {

        if (selectedTracks.isEmpty()) {
            return;
        }

        Color currentSelection = selectedTracks.iterator().next().getColor();

        Color color = UIUtilities.showColorChooserDialog(
                "Select Track Color (Negative Values)",
                currentSelection);

        if (color == null) {
            return;
        }

        for (Track track : selectedTracks) {
            track.setAltColor(ColorUtilities.modifyAlpha(color, currentSelection.getAlpha()));
        }
        refresh();

    }

    public static void exportTrackNames(final Collection<Track> selectedTracks) {

        if (selectedTracks.isEmpty()) {
            return;
        }

        File file = FileDialogUtils.chooseFile("Export track names",
                PreferenceManager.getInstance().getLastTrackDirectory(),
                new File("trackNames.tab"),
                FileDialogUtils.SAVE);

        if (file == null) {
            return;
        }

        PrintWriter pw = null;
        try {
             pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));

            List<String> attributes = AttributeManager.getInstance().getVisibleAttributes();

            pw.print("Name");
            for(String att : attributes) {
                pw.print("\t" + att);
            }
            pw.println();

            for (Track track : selectedTracks) {
                //We preserve the alpha value. This is motivated by MergedTracks
                pw.print(track.getName());

                for(String att : attributes) {
                    String val = track.getAttributeValue(att);
                    pw.print("\t" + (val == null ? "" : val));
                }
                pw.println();
            }


        } catch (IOException e) {
            MessageUtils.showErrorMessage("Error writing to file", e);
            log.error(e);
        } finally {
            if(pw != null) pw.close();
        }

    }


    public static JMenuItem getCopyDetailsItem(final Feature f, final TrackClickEvent evt) {
        JMenuItem item = new JMenuItem("Copy Details to Clipboard");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                ReferenceFrame frame = evt.getFrame();
                int mouseX = evt.getMouseEvent().getX();

                double location = frame.getChromosomePosition(mouseX);
                if (f instanceof IGVFeature) {
                    String details = ((IGVFeature) f).getValueString(location, null);
                    if (details != null) {
                        details = details.replace("<br>", System.getProperty("line.separator"));
                        details += System.getProperty("line.separator") +
                                f.getChr() + ":" + (f.getStart() + 1) + "-" + f.getEnd();
                        StringUtils.copyTextToClipboard(details);
                    }
                }
            }
        });
        return item;
    }


    public static JMenuItem getCopySequenceItem(final Feature f) {
        JMenuItem item = new JMenuItem("Copy Sequence");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                IGV.copySequenceToClipboard(genome, f.getChr(), f.getStart(), f.getEnd());
            }
        });
        return item;
    }


    public static JMenuItem getBlatItem(final Feature f) {
        JMenuItem item = new JMenuItem("Blat Sequence");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                BlatClient.doBlatQuery(f.getChr(), f.getStart(), f.getEnd());
            }
        });
        return item;
    }


    /**
     * Return a representative track height to use as the default.  For now
     * using the median track height.
     *
     * @return
     */
    public static int getRepresentativeTrackHeight(Collection<Track> tracks) {

        double[] heights = new double[tracks.size()];
        int i = 0;
        for (Track track : tracks) {
            heights[i] = track.getHeight();
            i++;
        }
        int medianTrackHeight = (int) Math.round(StatUtils.percentile(heights, 50));
        if (medianTrackHeight > 0) {
            return medianTrackHeight;
        }

        return PreferenceManager.getInstance().getAsInt(PreferenceManager.INITIAL_TRACK_HEIGHT);

    }

    public static void refresh() {
        if (IGV.hasInstance()) {
            IGV.getInstance().showLoadedTrackCount();
            IGV.getInstance().doRefresh();
        }
    }

    public static JMenuItem getChangeTrackHeightItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Change Track Height...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                changeTrackHeight(selectedTracks);
            }
        });
        return item;
    }

    public static JMenuItem getChangeKMPlotItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Kaplan-Meier Plot...");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {

                // If one or fewer tracks are selected assume the intent is to use all tracks.  A right-click
                // will always result in one selected track.
                Collection<Track> tracks = selectedTracks.size() > 1 ? selectedTracks :
                        IGV.getInstance().getAllTracks();
                KMPlotFrame frame = new KMPlotFrame(tracks);
                frame.setVisible(true);
            }
        });

        // The Kaplan-Meier plot requires sample information, specifically survival, sample, and censure.  We
        // can't know if these columns exist, but we can at least know if sample-info has been loaded.
        // 3-4 columns always exist by default, more indicate at least some sample attributes are defined.
        boolean sampleInfoLoaded = AttributeManager.getInstance().getAttributeNames().size() > 4;
        item.setEnabled(sampleInfoLoaded);
        return item;
    }

    public static JMenuItem getChangeFeatureWindow(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Set Feature Visibility Window...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                changeFeatureVisibilityWindow(selectedTracks);
            }
        });
        return item;
    }

    public static JMenuItem getChangeFontSizeItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Change Font Size...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                changeFontSize(selectedTracks);
            }
        });
        return item;
    }


    // Experimental methods follow

    public static JMenuItem getShowSortFramesItem(final Track track) {

        final JCheckBoxMenuItem item = new JCheckBoxMenuItem("Sort frames");

        item.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                Runnable runnable = new Runnable() {
                    public void run() {
                        FrameManager.sortFrames(track);
                        IGV.getInstance().resetFrames();
                    }
                };
                LongRunningTask.submit(runnable);
            }

        });
        return item;
    }

}

