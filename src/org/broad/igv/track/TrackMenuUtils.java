/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.track;


import com.jidesoft.swing.JidePopupMenu;
import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.stats.KMPlotTest;
import org.broad.tribble.Feature;
import org.broad.igv.feature.SequenceManager;
import org.broad.igv.renderer.*;
import org.broad.igv.ui.DataRangeDialog;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.HeatmapScaleDialog;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * @author jrobinso
 */
public class TrackMenuUtils {

    static Logger log = Logger.getLogger(TrackMenuUtils.class);
    final static String LEADING_HEADING_SPACER = "  ";
    static Font boldFont = FontManager.getScalableFont(Font.BOLD, 12);
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


    public static JidePopupMenu getEmptyPopup(String title) {
        JidePopupMenu menu = new JidePopupMenu();

        JLabel popupTitle = new JLabel(LEADING_HEADING_SPACER + title, JLabel.CENTER);
        popupTitle.setFont(boldFont);
        if (popupTitle != null) {
            menu.add(popupTitle);
            menu.addSeparator();
        }
        return menu;
    }

    /**
     * Return a popup menu with items applicable to the collection of tracks.
     *
     * @param tracks
     * @return
     */
    public static JPopupMenu getPopupMenu(Collection<Track> tracks, String title, TrackClickEvent te) {

        if (log.isDebugEnabled()) {
            log.debug("enter getPopupMenu");
        }

        JidePopupMenu menu = getEmptyPopup(title);
        addStandardItems(menu, tracks, te);
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
        if (dataTracksOnly) {
            menu.addSeparator();
            addDataItems(menu, tracks);
        } else if (featureTracksOnly) {
            menu.addSeparator();
            addFeatureItems(menu, tracks, te);
        }

        ReferenceFrame frame = te.getFrame();
        if (frame != null) {
            menu.addSeparator();
            addZoomItems(menu, frame);
        }

        menu.addSeparator();
        menu.add(getRemoveMenuItem(tracks));

    }

    public static void addZoomItems(JPopupMenu menu, final ReferenceFrame frame) {

        if (FrameManager.isGeneListMode()) {
            JMenuItem item = new JMenuItem("Reset view");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    frame.reset();
                    // TODO -- paint only panels for this frame
                }
            });
            menu.add(item);
        }

        JMenuItem zoomInItem = new JMenuItem("Zoom in");
        zoomInItem.addActionListener(new TrackActionListener() {

            public void action() {
                frame.incrementZoom(1);
            }
        });
        menu.add(zoomInItem);

        JMenuItem zoomOutItem = new JMenuItem("Zoom out");
        zoomOutItem.addActionListener(new TrackActionListener() {

            public void action() {
                frame.incrementZoom(-1);
            }
        });
        menu.add(zoomOutItem);


    }


    /**
     * Return popup menu with items applicable to data tracks
     *
     * @return
     */
    private static void addDataItems(JPopupMenu menu, final Collection<Track> tracks) {

        if (log.isDebugEnabled()) {
            log.debug("enter getDataPopupMenu");
        }

        final String[] labels = {"Heatmap", "Bar Chart", "Scatterplot", "Line Plot"};
        final Class[] renderers = {HeatmapRenderer.class, BarChartRenderer.class,
                ScatterplotRenderer.class, LineplotRenderer.class
        };

        //JLabel popupTitle = new JLabel(LEADING_HEADING_SPACER + title, JLabel.CENTER);

        JLabel rendererHeading = new JLabel(LEADING_HEADING_SPACER + "Type of Graph", JLabel.LEFT);
        rendererHeading.setFont(boldFont);

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
            item.addActionListener(new TrackActionListener() {
                public void action() {
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

        if (!avaibleWindowFunctions.isEmpty() || !currentWindowFunctions.isEmpty()) {
            JLabel statisticsHeading = new JLabel(LEADING_HEADING_SPACER + "Windowing Function", JLabel.LEFT);
            statisticsHeading.setFont(boldFont);

            menu.add(statisticsHeading);

            for (final WindowFunction wf : ORDERED_WINDOW_FUNCTIONS) {
                JCheckBoxMenuItem item = new JCheckBoxMenuItem(wf.getDisplayName());
                if (avaibleWindowFunctions.contains(wf) || currentWindowFunctions.contains(wf)) {
                    if (currentWindowFunctions.contains(wf)) {
                        item.setSelected(true);
                    }
                    item.addActionListener(new TrackActionListener() {

                        public void action() {
                            changeStatType(wf.toString(), tracks);
                        }
                    });
                    menu.add(item);
                }
            }
            menu.addSeparator();
        }

        JLabel scaleHeading = new JLabel(LEADING_HEADING_SPACER + "Data Range", JLabel.LEFT);
        scaleHeading.setFont(boldFont);
        menu.add(scaleHeading);

        menu.add(getDataRangeItem(tracks));
        menu.add(getHeatmapScaleItem(tracks));

        if (tracks.size() > 0) {
            menu.add(getLogScaleItem(tracks));
        }

        menu.add(getAutoscaleItem(tracks));


        menu.add(getShowDataRangeItem(tracks));


        //menu.add(getChangeKMPlotItem(tracks));

    }

    /**
     * Return popup menu with items applicable to feature tracks
     *
     * @return
     */
    private static void addFeatureItems(JPopupMenu featurePopupMenu, final Collection<Track> tracks, TrackClickEvent te) {


        featurePopupMenu.add(getExpandCollapseItem(tracks));

        featurePopupMenu.add(getChangeFontSizeItem(tracks));

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
    }

    /**
     * Popup menu with items applicable to both feature and data tracks
     *
     * @return
     */
    private static void addSharedItems(JPopupMenu menu, final Collection<Track> tracks, boolean hasFeatureTracks) {

        //JLabel trackSettingsHeading = new JLabel(LEADING_HEADING_SPACER + "Track Settings", JLabel.LEFT);
        //trackSettingsHeading.setFont(boldFont);
        //menu.add(trackSettingsHeading);

        menu.add(getTrackRenameItem(tracks));

        String colorLabel = hasFeatureTracks
                ? "Change Track Color" : "Change Track Color (Positive Values)";
        JMenuItem item = new JMenuItem(colorLabel);
        item.addActionListener(new TrackActionListener() {
            public void action() {
                changeTrackColor(tracks);
            }
        });
        menu.add(item);

        if (!hasFeatureTracks) {

            // Change track color by attribute
            item = new JMenuItem("Change Track Color (Negative Values)");
            item.setToolTipText(
                    "Change the alternate track color.  This color is used when graphing negative values");
            item.addActionListener(new TrackActionListener() {
                public void action() {
                    changeAltTrackColor(tracks);
                }
            });
            menu.add(item);
        }

        menu.add(getChangeTrackHeightItem(tracks));
    }


    private static void changeStatType(String statType, Collection<Track> selectedTracks) {
        for (Track track : selectedTracks) {
            track.setStatType(WindowFunction.valueOf(statType));
        }
        refresh();
    }


    public static JMenuItem getTrackRenameItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Rename Track");
        item.addActionListener(new TrackActionListener() {

            public void action() {
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

        item.addActionListener(new TrackActionListener() {

            public void action() {
                if (selectedTracks.size() > 0) {

                    ContinuousColorScale colorScale = selectedTracks.iterator().next().getColorScale();
                    HeatmapScaleDialog dlg = new HeatmapScaleDialog(IGVMainFrame.getInstance(), colorScale);

                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        colorScale = dlg.getColorScale();

                        // dlg.isFlipAxis());
                        for (Track track : selectedTracks) {
                            track.setColorScale(colorScale);
                        }
                        IGVMainFrame.getInstance().repaint();
                    }

                }

            }
        });
        return item;
    }

    public static JMenuItem getDataRangeItem(final Collection<Track> selectedTracks) {
        JMenuItem item = new JMenuItem("Set Data Range...");

        item.addActionListener(new TrackActionListener() {

            public void action() {
                if (selectedTracks.size() > 0) {
                    // Create a datarange that spans the extent of prev tracks range
                    float mid = 0;
                    float min = Float.MAX_VALUE;
                    float max = Float.MIN_VALUE;
                    boolean drawBaseline = false;
                    for (Track t : selectedTracks) {
                        DataRange dr = t.getDataRange();
                        min = Math.min(min, dr.getMinimum());
                        max = Math.max(max, dr.getMaximum());
                        mid += dr.getBaseline();
                    }
                    mid /= selectedTracks.size();
                    if (mid < min) {
                        mid = min;
                    } else if (mid > max) {
                        min = max;
                    }

                    DataRange prevAxisDefinition = new DataRange(min, mid, max, drawBaseline);
                    DataRangeDialog dlg = new DataRangeDialog(IGVMainFrame.getInstance(), prevAxisDefinition);
                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        min = Math.min(dlg.getMax(), dlg.getMin());
                        max = Math.max(dlg.getMin(), dlg.getMax());
                        mid = dlg.getBase();
                        mid = Math.max(min, Math.min(mid, max));

                        DataRange axisDefinition = new DataRange(dlg.getMin(), dlg.getBase(), dlg.getMax());

                        for (Track track : selectedTracks) {
                            track.setDataRange(axisDefinition);
                            if (track instanceof DataTrack) {
                                ((DataTrack) track).setAutoscale(false);
                            }
                        }
                        IGVMainFrame.getInstance().repaint();
                    }

                }
            }
        });
        return item;
    }

    private static JMenuItem getLogScaleItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute


        final JCheckBoxMenuItem logScaleItem = new JCheckBoxMenuItem("Log scale");
        final boolean logScale = selectedTracks.iterator().next().getDataRange().isLog();
        logScaleItem.setSelected(logScale);
        logScaleItem.addActionListener(new TrackActionListener() {

            public void action() {
                DataRange.Type scaleType = logScaleItem.isSelected() ?
                        DataRange.Type.LOG :
                        DataRange.Type.LINEAR;
                for (Track t : selectedTracks) {
                    t.getDataRange().setType(scaleType);
                }
                IGVMainFrame.getInstance().repaintDataPanels();
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
                if (t instanceof DataTrack && ((DataTrack) t).isAutoscale()) {
                    autoScale = true;
                    break;
                }
            }

            autoscaleItem.setSelected(autoScale);
            autoscaleItem.addActionListener(new TrackActionListener() {

                public void action() {

                    boolean autoScale = autoscaleItem.isSelected();
                    for (Track t : selectedTracks) {
                        if (t instanceof DataTrack) {
                            ((DataTrack) t).setAutoscale(autoScale);
                        }
                    }
                    IGVMainFrame.getInstance().repaintDataPanels();
                }
            });
        }
        return autoscaleItem;
    }

    private static JMenuItem getShowDataRangeItem(final Collection<Track> selectedTracks) {

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
            item.addActionListener(new TrackActionListener() {

                public void action() {

                    boolean showDataRange = item.isSelected();
                    for (Track t : selectedTracks) {
                        if (t instanceof DataTrack) {
                            ((DataTrack) t).setShowDataRange(showDataRange);
                        }
                    }
                    IGVMainFrame.getInstance().repaintDataPanels();
                }
            });
        }
        return item;
    }


    public static JMenuItem getExpandCollapseItem(final Collection<Track> tracks) {

        // If any tracks are expanded show the "Collapse" option, otherwise expand
        boolean expanded = false;
        for (Track track : tracks) {
            if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
                expanded = true;
                break;
            }
        }

        // Show Multi-Level Features
        String text = expanded ? "Collapse Track" : "Expand Track";
        if (tracks.size() > 1) {
            text += "s";
        }
        JMenuItem item = new JMenuItem(text);
        item.addActionListener(new TrackActionListener() {

            public void action() {
                toggleExpandedState(tracks);
                IGVMainFrame.getInstance().doRefresh();
            }
        });

        return item;
    }

    public static JMenuItem getRemoveMenuItem(final Collection<Track> selectedTracks) {

        boolean multiple = selectedTracks.size() > 1;

        JMenuItem item = new JMenuItem("Remove Track" + (multiple ? "s" : ""));
        item.addActionListener(new TrackActionListener() {

            public void action() {
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
                JDialog dialog = optionPane.createDialog(IGVMainFrame.getInstance(), "Remove The Following Tracks");
                dialog.setVisible(true);

                Object choice = optionPane.getValue();
                if ((choice == null) || (JOptionPane.YES_OPTION != ((Integer) choice).intValue())) {
                    return;
                }

                IGVMainFrame.getInstance().getTrackManager().removeTracks(selectedTracks);
                IGVMainFrame.getInstance().updateTrackState();
            }
        });
        return item;
    }

    /**
     * Toggle selected multi-level tracks
     * if (track instanceof FeatureTrack) {
     * if (((FeatureTrack) track).getDisplayMode() == Track.DisplayMode.EXPANDED) {
     * expanded = true;
     * break;
     * }
     * }
     */
    public static void toggleExpandedState(final Collection<Track> selectedTracks) {

        for (Track track : selectedTracks) {
            boolean expanded = track.getDisplayMode() == Track.DisplayMode.EXPANDED;
            if (!expanded) {
                track.setDisplayMode(Track.DisplayMode.EXPANDED);
            } else {
                track.setDisplayMode(Track.DisplayMode.COLLAPSED);
            }
        }
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
        String newName = JOptionPane.showInputDialog(IGVMainFrame.getInstance(), "Enter new name: ", t.getName());

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
        int origValueKB = Math.max(1, origValue / 1000);
        int value = getIntValue("Visibility window (kb)", origValueKB);
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


    private static int getIntValue(String parameter, int value) {

        while (true) {

            String height = JOptionPane.showInputDialog(
                    IGVMainFrame.getInstance(), parameter + ": ",
                    String.valueOf(value));

            if ((height == null) || height.trim().equals("")) {
                return Integer.MIN_VALUE;   // <= the logical "null" value
            }

            try {
                value = Integer.parseInt(height);
                return value;
            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGVMainFrame.getInstance(),
                        parameter + " must be an integer number.");
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
        item.addActionListener(new TrackActionListener() {

            public void action() {

                ReferenceFrame frame = evt.getFrame();
                int mouseX = evt.getMouseEvent().getX();

                double location = frame.getChromosomePosition(mouseX);
                if (f instanceof IGVFeature) {
                    String details = ((IGVFeature) f).getValueString(location, null);
                    if (details != null) {
                        details = details.replace("<br>", System.getProperty("line.separator"));
                        details += System.getProperty("line.separator") +
                                f.getChr() + ":" + (f.getStart() + 1) + "-" + f.getEnd();
                        StringSelection stringSelection = new StringSelection(details);
                        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                        clipboard.setContents(stringSelection, null);
                    }
                }
            }
        });
        return item;
    }


    public static JMenuItem getCopySequenceItem(final Feature f) {
        JMenuItem item = new JMenuItem("Copy Sequence");
        item.addActionListener(new TrackActionListener() {

            public void action() {
                String genomeId = GenomeManager.getInstance().getGenomeId();
                String chr = f.getChr();
                int start = f.getStart();
                int end = f.getEnd();
                byte[] seqBytes = SequenceManager.readSequence(genomeId, chr, start, end);
                if (seqBytes == null) {
                    MessageUtils.showMessage("Sequence not available");
                } else {
                    String sequence = new String(seqBytes);
                    StringSelection stringSelection = new StringSelection(sequence);
                    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                    clipboard.setContents(stringSelection, null);
                }
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
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {

                IGVMainFrame.getInstance().showLoadedTrackCount();
                IGVMainFrame.getInstance().getContentPane().repaint();
            }
        });
    }

    public static JMenuItem getChangeTrackHeightItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Change Track Height...");
        item.addActionListener(new TrackActionListener() {

            public void action() {
                changeTrackHeight(selectedTracks);
            }
        });
        return item;
    }

    public static JMenuItem getChangeKMPlotItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Kaplan-Meier Plot...");
        item.addActionListener(new TrackActionListener() {

            public void action() {
                KMPlotTest.openPlot(selectedTracks);
            }
        });
        return item;
    }

    public static JMenuItem getChangeFeatureWindow(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Set Feature Visibility Window...");
        item.addActionListener(new TrackActionListener() {

            public void action() {
                changeFeatureVisibilityWindow(selectedTracks);
            }
        });
        return item;
    }

    public static JMenuItem getChangeFontSizeItem(final Collection<Track> selectedTracks) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Change Font Size");
        item.addActionListener(new TrackActionListener() {

            public void action() {
                changeFontSize(selectedTracks);
            }
        });
        return item;
    }

    public static abstract class TrackActionListener implements ActionListener {

        public abstract void action();

        public void actionPerformed(ActionEvent actionEvent) {
            action();
            IGVMainFrame.getInstance().getTrackManager().clearSelections();

        }
    }
}

