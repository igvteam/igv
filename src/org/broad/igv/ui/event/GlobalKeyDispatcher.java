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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.event;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.charts.ScatterPlotUtils;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.variant.VariantTrack;
import htsjdk.tribble.Feature;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public class GlobalKeyDispatcher implements KeyEventDispatcher {

    private static Logger log = Logger.getLogger(GlobalKeyDispatcher.class);

    private final InputMap inputMap = new InputMap();
    private final ActionMap actionMap = new ActionMap();

    public GlobalKeyDispatcher() {
        init();
    }

    public boolean dispatchKeyEvent(KeyEvent event) {

        KeyStroke ks = KeyStroke.getKeyStrokeForEvent(event);
        String actionKey = (String) inputMap.get(ks);

        // Disable tooltip if any modifier control key is pressed
        if (event.getKeyCode() == KeyEvent.VK_CONTROL || event.getKeyCode() == KeyEvent.VK_ALT) {
            boolean flag = !(event.isControlDown() || event.isAltDown() || event.isMetaDown());
            ToolTipManager.sharedInstance().setEnabled(flag);
        }

        if (event.getKeyCode() == KeyEvent.VK_ESCAPE) {
            IGV.getInstance().clearSelections();
            IGV.getInstance().repaint();
            return true;
        }


        if (actionKey != null) {
            Action action = actionMap.get(actionKey);
            if (action != null && action.isEnabled()) {
                // I'm not sure about the parameters
                action.actionPerformed(
                        new ActionEvent(event.getSource(), event.getID(),
                                actionKey, ((KeyEvent) event).getModifiers()));
                return true; // consume event
            }
        }

        return false;
    }

    /**
     * Initialize the input and action map.   The indirection here strikes me as odd but it is apparently the standard pattern.
     */
    public void init() {

        final IGV igv = IGV.getInstance();
        final PreferenceManager prefMgr = PreferenceManager.getInstance();

        // Next feature
        final KeyStroke nextKey = KeyStroke.getKeyStroke(KeyEvent.VK_F, KeyEvent.CTRL_MASK, false);
        final Action nextAction = new EnableWrappedAction(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                nextFeature(true);
            }
        });
        inputMap.put(nextKey, "nextFeature");
        actionMap.put("nextFeature", nextAction);

        // Previous feature
        final KeyStroke prevKey = KeyStroke.getKeyStroke(KeyEvent.VK_B, KeyEvent.CTRL_MASK, false);
        final Action prevAction = new EnableWrappedAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                nextFeature(false);
            }
        });
        inputMap.put(prevKey, "prevFeature");
        actionMap.put("prevFeature", prevAction);

        // Next exon
        final KeyStroke nextExonKey = KeyStroke.getKeyStroke(KeyEvent.VK_F, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);
        final Action nextExonAction = new EnableWrappedAction(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                nextExon(true);
            }
        });
        inputMap.put(nextExonKey, "nextExon");
        actionMap.put("nextExon", nextExonAction);

        // Previous exon
        final KeyStroke prevExonKey = KeyStroke.getKeyStroke(KeyEvent.VK_B, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);
        final Action prevExonAction = new EnableWrappedAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                nextExon(false);
            }
        });
        inputMap.put(prevExonKey, "prevExon");
        actionMap.put("prevExon", prevExonAction);

        // Show extras menu
        final KeyStroke extrasKey = KeyStroke.getKeyStroke(KeyEvent.VK_T, KeyEvent.ALT_MASK, false);
        final Action extrasAction = new EnableWrappedAction(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                igv.enableExtrasMenu();
            }
        });
        inputMap.put(extrasKey, "tools");
        actionMap.put("tools", extrasAction);

        // Create region-on-interest
        final KeyStroke regionKey = KeyStroke.getKeyStroke(KeyEvent.VK_R, KeyEvent.CTRL_MASK, false);
        final Action regionAction = new EnableWrappedAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (FrameManager.isGeneListMode()) {
                    return;
                }
                Range currentRange = FrameManager.getDefaultFrame().getCurrentRange();
                RegionOfInterest regionOfInterest =
                        new RegionOfInterest(
                                currentRange.getChr(),
                                currentRange.getStart(),
                                currentRange.getEnd(),
                                null);
                igv.addRegionOfInterest(regionOfInterest);
            }
        });
        inputMap.put(regionKey, "region");
        actionMap.put("region", regionAction);

        // Create region-of-interest at center of view (1 bp wide)
        final KeyStroke regionCenterKey = KeyStroke.getKeyStroke(KeyEvent.VK_R, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);
        final Action regionCenterAction = new EnableWrappedAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (FrameManager.isGeneListMode()) {
                    return;
                }
                int center = (int) FrameManager.getDefaultFrame().getCenter();
                RegionOfInterest regionOfInterest =
                        new RegionOfInterest(
                                FrameManager.getDefaultFrame().getChrName(),
                                center,
                                center + 1,
                                null);
                igv.addRegionOfInterest(regionOfInterest);
            }
        });
        inputMap.put(regionCenterKey, "regionCenter");
        actionMap.put("regionCenter", regionCenterAction);

        // Sort alignments
        final KeyStroke sortByLastKey = KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK, false);
        final Action sorAlignmentTracksAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String sortOptionString = prefMgr.get(PreferenceManager.SAM_SORT_OPTION);
                if (sortOptionString != null) {
                    try {
                        AlignmentTrack.SortOption option = AlignmentTrack.SortOption.valueOf(sortOptionString);
                        String lastSortTag = prefMgr.get(PreferenceManager.SAM_SORT_BY_TAG);

                        igv.sortAlignmentTracks(option, lastSortTag);
                        igv.repaintDataPanels();
                    } catch (IllegalArgumentException e1) {
                        log.error("Unrecognized sort option: " + sortOptionString);
                    }
                }
            }
        };
        inputMap.put(sortByLastKey, "sortByLast");
        actionMap.put("sortByLast", sorAlignmentTracksAction);

        // Open scatter plot
        final KeyStroke scatterplotKey = KeyStroke.getKeyStroke(KeyEvent.VK_P, KeyEvent.CTRL_MASK, false);
        final Action scatterplotAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (ScatterPlotUtils.hasPlottableTracks()) {
                    ReferenceFrame defaultFrame = FrameManager.getDefaultFrame();
                    String chr = defaultFrame.getChrName();
                    int start = (int) defaultFrame.getOrigin();
                    int end = (int) defaultFrame.getEnd();
                    int zoom = defaultFrame.getZoom();
                    ScatterPlotUtils.openPlot(chr, start, end, zoom);
                }
            }
        };
        inputMap.put(scatterplotKey, "scatterPlot");
        actionMap.put("scatterPlot", scatterplotAction);

        // Back button
        final KeyStroke backKey1 = KeyStroke.getKeyStroke(KeyEvent.VK_CLOSE_BRACKET, KeyEvent.META_DOWN_MASK, false);
        final KeyStroke backKey2 = KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, KeyEvent.ALT_DOWN_MASK, false);
        final Action backAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                igv.getSession().getHistory().back();
            }
        };
        inputMap.put(backKey1, "back");

        inputMap.put(backKey2, "back");
        actionMap.put("back", backAction);

        // Forward button
        final KeyStroke forwardKey1 = KeyStroke.getKeyStroke(KeyEvent.VK_OPEN_BRACKET, KeyEvent.META_DOWN_MASK, false);
        final KeyStroke forwardKey2 = KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, KeyEvent.ALT_DOWN_MASK, false);
        final Action forwardAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                igv.getSession().getHistory().forward();
            }
        };
        inputMap.put(forwardKey1, "forward");
        inputMap.put(forwardKey2, "forward");

        // Toggle alignment "complete read only" option
        final KeyStroke completeReadKey = KeyStroke.getKeyStroke(KeyEvent.VK_1, KeyEvent.CTRL_DOWN_MASK, false);
        final Action completeReadAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                boolean currentSetting = prefMgr.getAsBoolean(PreferenceManager.SAM_COMPLETE_READS_ONLY);
                prefMgr.put(PreferenceManager.SAM_COMPLETE_READS_ONLY, !currentSetting);
                igv.repaintDataPanels();
            }
        };
        inputMap.put(completeReadKey, "completeReads");
        actionMap.put("completeReads", completeReadAction);


    }

    /**
     * Move to the next exon in the feature located at the center, if:
     * -there is such a feature
     * -the track is expanded
     * -a single feature row is selected
     * -the feature has multiple exons
     * -there is an exon forward or backward to jump to
     *
     * @param forward
     */
    private void nextExon(boolean forward) {

        // Ignore (Disable) if we are in gene list mode
        if (FrameManager.isGeneListMode()) {
            return;
        }

        ReferenceFrame vc = FrameManager.getDefaultFrame();
        Collection<Track> tracks = IGV.getInstance().getSelectedTracks();
        if (tracks.size() == 1) {
            Track t = tracks.iterator().next();
            if (!(t instanceof FeatureTrack)) {
                //JOptionPane.showMessageDialog(IGV.getInstance(),
                //        "Track panning is not enabled for data tracks.");
                return;
            }

            Exon e = null;
            if (t instanceof FeatureTrack) {
                int center = (int) vc.getCenter();
                FeatureTrack ft = (FeatureTrack) t;
                if (ft.getDisplayMode() == Track.DisplayMode.COLLAPSED ||
                        ft.getSelectedFeatureRowIndex() == FeatureTrack.NO_FEATURE_ROW_SELECTED) {
                    MessageUtils.showMessage(
                            "Exon navigation is only allowed when track is expanded and a single " +
                                    "feature row is selected.");
                    return;
                }
                List<Feature> featureList = ft.getFeaturesAtPositionInFeatureRow(center, ft.getSelectedFeatureRowIndex(), vc);
                Feature feature = featureList != null && featureList.size() > 0 ? featureList.get(0) : null;

                if (feature == null)
                    return;
                if (feature instanceof BasicFeature) {
                    BasicFeature bf = (BasicFeature) feature;
                    java.util.List<Exon> exons = bf.getExons();
                    if (exons == null || exons.isEmpty()) {
                        MessageUtils.showMessage("At least one centered feature does not have exon structure");
                        return;
                    }

                    if (forward) {
                        for (Exon exon : bf.getExons()) {
                            //the "+ 1" here is necessary because the rounding in the recentering method
                            //sometimes places the center one base off.  This should be perfectly safe,
                            //but it does assume no one's abusing the exon datastructure and creating
                            //exons that are right next to each other.
                            if (exon.getStart() > vc.getCenter() + 1) {
                                e = exon;
                                break;
                            }
                        }
                    } else {
                        for (int i = exons.size() - 1; i >= 0; i--) {
                            Exon exon = exons.get(i);
                            if (exon.getEnd() < vc.getCenter()) {
                                e = exon;
                                break;
                            }
                        }
                    }
                }

                if (e != null) {
                    vc.centerOnLocation(forward ? e.getStart() : e.getEnd());
                    int i = 3;
                }

            }
            //todo: implement handling for VariantTrack


        } else {
            MessageUtils.showMessage("To use track panning you must first select a single feature track.");
        }
    }


    /**
     * Skip to the next feature in the selected track.
     *
     * @param forward the direction, true for forward and false for back
     */
    private void nextFeature(boolean forward) {

        // Ignore (Disable) if we are in gene list mode
        if (FrameManager.isGeneListMode()) {
            return;
        }

        ReferenceFrame frame = FrameManager.getDefaultFrame();
        Collection<Track> tracks = IGV.getInstance().getSelectedTracks();
        if (tracks.size() == 1) {
            try {
                Track t = tracks.iterator().next();
                if (!(t instanceof FeatureTrack || t instanceof VariantTrack)) {
                    //JOptionPane.showMessageDialog(IGV.getInstance(),
                    //        "Track panning is not enabled for data tracks.");
                    return;
                }

                Feature f = null;
                if (t instanceof FeatureTrack) {
                    f = ((FeatureTrack) t).nextFeature(frame.getChrName(), frame.getCenter(), forward, frame);
                } else if (t instanceof VariantTrack) {
                    f = ((VariantTrack) t).nextFeature(frame.getChrName(), frame.getCenter(), forward, frame);
                }

                if (f != null) {
                    String chr = GenomeManager.getInstance().getCurrentGenome().getChromosomeAlias(f.getChr());
                    double newCenter = f.getStart();
                    if (!chr.equals(frame.getChrName())) {
                        // Switch chromosomes.  We have to do some tricks to maintain the same resolution scale.
                        double range = frame.getEnd() - frame.getOrigin();
                        int newOrigin = (int) Math.max(newCenter - range / 2, 0);
                        int newEnd = (int) (newOrigin + range);
                        frame.jumpTo(chr, newOrigin, newEnd);
                    } else {
                        frame.centerOnLocation(newCenter);
                    }
                }
            } catch (IOException e) {
                MessageUtils.showErrorMessage("Error encountered reading features: " + e.getMessage(), e);

            }
        } else {
            MessageUtils.showMessage("To use track panning you must first select a single feature track.");
        }


    }

    /**
     * TODO I'm actually pretty sure this class doesn't do what it's intended to do,
     * but I just refactored it to condense code, there were no functional changes.
     * -JS
     */
    private class EnableWrappedAction extends AbstractAction {

        private Action action;

        private EnableWrappedAction(Action action) {
            this.action = action;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            setEnabled(false); // stop any other events from interfering
            this.action.actionPerformed(e);
            setEnabled(true);

        }
    }

}
