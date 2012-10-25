/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;


import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.*;
import org.broad.igv.ui.DataRangeDialog;
import org.broad.igv.ui.HeatmapScaleDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.stats.KMPlotFrame;
import org.broad.tribble.Feature;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;

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

        //JidePopupMenu menu = new JidePopupMenu();
        IGVPopupMenu menu = new IGVPopupMenu();

        JLabel popupTitle = new JLabel(LEADING_HEADING_SPACER + title, JLabel.CENTER);
        popupTitle.setFont(UIConstants.boldFont);
        if (popupTitle != null) {
            menu.add(popupTitle);
            menu.addSeparator();
        }

        addStandardItems(menu, tracks, te);

        menu.addPopupMenuListener(new PopupMenuListener() {
            public void popupMenuWillBecomeVisible(PopupMenuEvent popupMenuEvent) {

            }

            public void popupMenuWillBecomeInvisible(PopupMenuEvent popupMenuEvent) {
                close();
            }

            public void popupMenuCanceled(PopupMenuEvent popupMenuEvent) {
                close();
            }

            private void close() {
                if (IGV.hasInstance()) {
                    IGV.getInstance().clearSelections();
                    IGV.getInstance().repaint();
                }
            }

        });

        return menu;

    }

    public static void addStandardItems(JPopupMenu menu, Collection<Track> tracks, TrackClickEvent te) {

        boolean hasDataTracks = false;
        boolean hasFeatureTracks = false;
        boolean hasOtherTracks = false;
        for (Track track : tracks) {

            //  TODO -- this is ugly, refactor to remove instanceof
            if (track instanceof DataTrack) {
                hasDataTracks = true;
            } else if (track instanceof FeatureTrack) {
                hasFeatureTracks = true;
            } else {
                hasOtherTracks = true;
            }
            if (hasDataTracks && hasFeatureTracks && hasOtherTracks) {
                break;
            }
        }

        boolean featureTracksOnly = hasFeatureTracks && !hasDataTracks && !hasOtherTracks;
        boolean dataTracksOnly = !hasFeatureTracks && hasDataTracks && !hasOtherTracks;

        addSharedItems(menu, tracks, hasFeatureTracks);
        menu.addSeparator();
        if (dataTracksOnly) {
            addDataItems(menu, tracks);
        } else if (featureTracksOnly) {
            addFeatureItems(menu, tracks, te);
        }

        menu.addSeparator();
        menu.add(getRemoveMenuItem(tracks));

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
                frame.incrementZoom(-1);
            }
        });
        menu.add(zoomOutItem);

        JMenuItem zoomInItem = new JMenuItem("Zoom in");
        zoomInItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                frame.incrementZoom(1);
            }
        });
        menu.add(zoomInItem);


    }


    /**
     * Return popup menu with items applicable to data tracks
     *
     * @return
     */
    public static void addDataItems(JPopupMenu menu, final Collection<Track> tracks) {

        if (log.isDebugEnabled()) {
            log.debug("enter getDataPopupMenu");
        }

        final String[] labels = {"Heatmap", "Bar Chart", "Points", "Line Plot"};
        final Class[] renderers = {HeatmapRenderer.class, BarChartRenderer.class,
                PointsRenderer.class, LineplotRenderer.class
        };

        //JLabel popupTitle = new JLabel(LEADING_HEADING_SPACER + title, JLabel.CENTER);

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

        // Create and renderer menu items
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


        menu.add(getDataRangeItem(tracks));
        menu.add(getHeatmapScaleItem(tracks));

        if (tracks.size() > 0) {
            menu.add(getLogScaleItem(tracks));
        }

        menu.add(getAutoscaleItem(tracks));


        menu.add(getShowDataRangeItem(tracks));

        menu.addSeparator();
        menu.add(getChangeKMPlotItem(tracks));

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
            }
        }

        featurePopupMenu.addSeparator();
        featurePopupMenu.add(getChangeFeatureWindow(tracks));

        //---------------------//
        //Track analysis
        if (Globals.toolsMenuEnabled && tracks.size() >= 2) {

            JMenuItem item = new JMenuItem("Create Overlap Track");
            item.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    String newName = JOptionPane.showInputDialog(IGV.getMainFrame(), "Enter overlap track name: ", "Overlaps");

                    if (newName == null || newName.trim() == "") {
                        return;
                    }
                    CombinedFeatureSource source = new CombinedFeatureSource(tracks, CombinedFeatureSource.Operation.MULTIINTER);
                    Track newTrack = new FeatureTrack("", newName, source);

                    IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTrack(newTrack);
                    IGV.getInstance().repaint();
                }
            });
            item.setEnabled(CombinedFeatureSource.checkBEDToolsPathValid());
            featurePopupMenu.add(item);
        }

        //--------------------//
    }


    /**
     * Popup menu with items applicable to both feature and data tracks
     *
     * @return
     */
    public static void addSharedItems(JPopupMenu menu, final Collection<Track> tracks, boolean hasFeatureTracks) {

        //JLabel trackSettingsHeading = new JLabel(LEADING_HEADING_SPACER + "Track Settings", JLabel.LEFT);
        //trackSettingsHeading.setFont(boldFont);
        //menu.add(trackSettingsHeading);

        menu.add(getTrackRenameItem(tracks));

        String colorLabel = hasFeatureTracks
                ? "Change Track Color..." : "Change Track Color (Positive Values)...";
        JMenuItem item = new JMenuItem(colorLabel);
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                changeTrackColor(tracks);
            }
        });
        menu.add(item);

        if (!hasFeatureTracks) {

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
                    float mid = 0;
                    float min = Float.MAX_VALUE;
                    float max = Float.MIN_VALUE;
                    boolean drawBaseline = true;
                    boolean isLog = true;
                    for (Track t : selectedTracks) {
                        DataRange dr = t.getDataRange();
                        min = Math.min(min, dr.getMinimum());
                        max = Math.max(max, dr.getMaximum());
                        mid += dr.getBaseline();
                        drawBaseline &= dr.isDrawBaseline();
                        isLog &= dr.isLog();
                    }
                    mid /= selectedTracks.size();
                    if (mid < min) {
                        mid = min;
                    } else if (mid > max) {
                        min = max;
                    }

                    DataRange prevAxisDefinition = new DataRange(min, mid, max, drawBaseline, isLog);
                    DataRangeDialog dlg = new DataRangeDialog(IGV.getMainFrame(), prevAxisDefinition);
                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        min = Math.min(dlg.getMax(), dlg.getMin());
                        max = Math.max(dlg.getMin(), dlg.getMax());
                        mid = dlg.getBase();
                        mid = Math.max(min, Math.min(mid, max));

                        DataRange axisDefinition = new DataRange(dlg.getMin(), mid, dlg.getMax(),
                                drawBaseline, dlg.isLog());

                        for (Track track : selectedTracks) {
                            track.setDataRange(axisDefinition);
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
            boolean autoScale = false;
            for (Track t : selectedTracks) {
                if (t instanceof DataTrack && ((DataTrack) t).isAutoScale()) {
                    autoScale = true;
                    break;
                }
            }

            autoscaleItem.setSelected(autoScale);
            autoscaleItem.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent evt) {

                    boolean autoScale = autoscaleItem.isSelected();
                    for (Track t : selectedTracks) {
                        if (t instanceof DataTrack) {
                            ((DataTrack) t).setAutoScale(autoScale);
                        }
                    }
                    IGV.getInstance().repaintDataPanels();
                }
            });
        }
        return autoscaleItem;
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
        boolean showAS = Boolean.parseBoolean(System.getProperty("showAS", "false"));
        if (showAS) {
            modes.put("Alternative Splice", Track.DisplayMode.ALTERNATIVE_SPLICE);
        }
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

                IGV.getInstance().removeTracks(selectedTracks);
                IGV.getInstance().doRefresh();
            }
        });
        return item;
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
        int value = getIntValue(parameter, getRepresentativeTrackHeight(selectedTracks));
        if (value == Integer.MIN_VALUE) {
            return;
        }

        value = Math.max(0, value);

        for (Track track : selectedTracks) {

            track.setHeight(value);

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
        double value = getNumericValue("Visibility window (kb)", origValueKB);
        if (value == Integer.MIN_VALUE) {
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
        int value = getIntValue(parameter, defaultValue);
        if (value == Integer.MIN_VALUE) {
            return;
        }

        for (Track track : selectedTracks) {
            track.setFontSize(value);
        }

        refresh();
    }


    public static int getIntValue(String parameter, int value) {

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

    public static double getNumericValue(String parameter, double value) {

        while (true) {

            String height = JOptionPane.showInputDialog(
                    IGV.getMainFrame(), parameter + ": ",
                    String.valueOf(value));

            if ((height == null) || height.trim().equals("")) {
                return Double.MIN_VALUE;   // <= the logical "null" value
            }

            try {
                value = Double.parseDouble(height);
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

        if (color != null) {
            for (Track track : selectedTracks) {
                track.setColor(color);
            }
            refresh();
        }

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

            track.setAltColor(color);
        }
        refresh();

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

}

