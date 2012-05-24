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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.event;

import org.broad.igv.charts.ScatterPlot;
import org.broad.igv.charts.ScatterPlotUtils;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.variant.VariantTrack;
import org.broad.tribble.Feature;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.IOException;
import java.util.Collection;

/**
 * @author jrobinso
 */
public class GlobalKeyDispatcher implements KeyEventDispatcher {

    private final InputMap keyStrokes = new InputMap();
    private final ActionMap actions = new ActionMap();

    public GlobalKeyDispatcher() {
        init();
    }

    public InputMap getInputMap() {
        return keyStrokes;
    }

    public ActionMap getActionMap() {
        return actions;
    }

    public boolean dispatchKeyEvent(KeyEvent event) {

        if (event.getKeyCode() == KeyEvent.VK_ESCAPE) {
            IGV.getInstance().clearSelections();
            IGV.getInstance().repaint();
            return true;
        }

        KeyStroke ks = KeyStroke.getKeyStrokeForEvent(event);
        String actionKey = (String) keyStrokes.get(ks);


        // Disable tooltip if any modifier control key is pressed
        if (event.getKeyCode() == KeyEvent.VK_CONTROL || event.getKeyCode() == KeyEvent.VK_ALT) {
            boolean flag = !(event.isControlDown() || event.isAltDown() || event.isMetaDown());
            ToolTipManager.sharedInstance().setEnabled(flag);
        }

        if (actionKey != null) {

            Action action = actions.get(actionKey);
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

    // Here for convenience,  move out of this class eventually

    public void init() {

        final KeyStroke nextKey = KeyStroke.getKeyStroke(KeyEvent.VK_F, KeyEvent.CTRL_MASK, false);
        final KeyStroke prevKey = KeyStroke.getKeyStroke(KeyEvent.VK_B, KeyEvent.CTRL_MASK, false);
        final KeyStroke toolsKey = KeyStroke.getKeyStroke(KeyEvent.VK_T, KeyEvent.ALT_MASK, false);
        final KeyStroke regionKey = KeyStroke.getKeyStroke(KeyEvent.VK_R, KeyEvent.CTRL_MASK, false);
        final KeyStroke regionCenterKey = KeyStroke.getKeyStroke(KeyEvent.VK_R, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);

        //dhmay adding 20101222
        final KeyStroke nextExonKey = KeyStroke.getKeyStroke(KeyEvent.VK_F, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);
        final KeyStroke prevExonKey = KeyStroke.getKeyStroke(KeyEvent.VK_B, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);


        final KeyStroke backKey1 = KeyStroke.getKeyStroke(KeyEvent.VK_CLOSE_BRACKET, KeyEvent.META_DOWN_MASK, false);
        final KeyStroke backKey2 = KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, KeyEvent.ALT_DOWN_MASK, false);
        final KeyStroke forwardKey1 = KeyStroke.getKeyStroke(KeyEvent.VK_OPEN_BRACKET, KeyEvent.META_DOWN_MASK, false);
        final KeyStroke forwardKey2 = KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, KeyEvent.ALT_DOWN_MASK, false);

        final KeyStroke statusWindowKey = KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_MASK, false);
        final KeyStroke scatterplotKey = KeyStroke.getKeyStroke(KeyEvent.VK_P, KeyEvent.CTRL_MASK, false);

        final Action toolAction = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                IGV.getInstance().enableExtrasMenu();
                setEnabled(true);
            }
        };

        final Action statusWindowAction = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                IGV.getInstance().openStatusWindow();
                setEnabled(true);
            }
        };

        final Action nextAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                nextFeature(true);
                setEnabled(true);
            }
        };
        final Action prevAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                nextFeature(false);
                setEnabled(true);
            }
        };

        //dhmay adding
        final Action nextExonAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                nextExon(true);
                setEnabled(true);
            }
        };
        final Action prevExonAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                nextExon(false);
                setEnabled(true);
            }
        };

        final Action regionAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                if (FrameManager.isGeneListMode()) {
                    return;
                }
                ReferenceFrame.Range currentRange = FrameManager.getDefaultFrame().getCurrentRange();
                RegionOfInterest regionOfInterest =
                        new RegionOfInterest(
                                currentRange.getChr(),
                                currentRange.getStart(),
                                currentRange.getEnd(),
                                null);
                // TODO -- get this ugly reference to IGV out of here
                IGV.getInstance().addRegionOfInterest(regionOfInterest);
                setEnabled(true);
            }
        };

        final Action regionCenterAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                setEnabled(false); // stop any other events from interfering
                if (FrameManager.isGeneListMode()) {
                    return;
                }
                int center = (int)  FrameManager.getDefaultFrame().getCenter();
                RegionOfInterest regionOfInterest =
                        new RegionOfInterest(
                                FrameManager.getDefaultFrame().getChrName(),
                                center,
                                center + 1,
                                null);
                // TODO -- get this ugly reference to IGV out of here
                IGV.getInstance().addRegionOfInterest(regionOfInterest);
                setEnabled(true);
            }
        };

        final Action backAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                IGV.getInstance().getSession().getHistory().back();
            }
        };

        final Action forwardAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                IGV.getInstance().getSession().getHistory().forward();
            }
        };

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

        getInputMap().put(nextKey, "nextFeature");
        getActionMap().put("nextFeature", nextAction);
        getInputMap().put(prevKey, "prevFeature");
        getActionMap().put("prevFeature", prevAction);

        //dhmay adding 20101222
        getInputMap().put(nextExonKey, "nextExon");
        getActionMap().put("nextExon", nextExonAction);
        getInputMap().put(prevExonKey, "prevExon");
        getActionMap().put("prevExon", prevExonAction);

        getInputMap().put(toolsKey, "tools");
        getActionMap().put("tools", toolAction);
        getInputMap().put(regionKey, "region");
        getActionMap().put("region", regionAction);
        getInputMap().put(regionCenterKey, "regionCenter");
        getActionMap().put("regionCenter", regionCenterAction);
        getInputMap().put(statusWindowKey, "statusWindow");
        getActionMap().put("statusWindow", statusWindowAction);
        getInputMap().put(scatterplotKey, "statusWindow");
        getActionMap().put("statusWindow", scatterplotAction);

        getInputMap().put(backKey1, "back");
        getInputMap().put(backKey2, "back");
        getActionMap().put("back", backAction);
        getInputMap().put(forwardKey1, "forward");
        getInputMap().put(forwardKey2, "forward");
        getActionMap().put("forward", forwardAction);

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
                Feature feature = ft.getFeatureAtPositionInFeatureRow(center, ft.getSelectedFeatureRowIndex(), vc);

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

        ReferenceFrame vc = FrameManager.getDefaultFrame();
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
                ReferenceFrame frame = FrameManager.getDefaultFrame();
                if (t instanceof FeatureTrack) {
                    f = ((FeatureTrack) t).nextFeature(vc.getChrName(), vc.getCenter(), forward, frame);
                } else if (t instanceof VariantTrack) {
                    f = ((VariantTrack) t).nextFeature(vc.getChrName(), vc.getCenter(), forward, frame);
                }

                if (f != null) {
                    String chr = GenomeManager.getInstance().getCurrentGenome().getChromosomeAlias(f.getChr());
                    double newCenter = f.getStart();
                    if (!chr.equals(vc.getChrName())) {
                        // Switch chromosomes.  We have to do some tricks to maintain the same resolution scale.
                        double range = vc.getEnd() - vc.getOrigin();
                        int newOrigin = (int) Math.max(newCenter - range / 2, 0);
                        int newEnd = (int) (newOrigin + range);
                        vc.jumpTo(chr, newOrigin, newEnd);
                    } else {
                        vc.centerOnLocation(newCenter);
                    }
                }
            } catch (IOException e) {
                //logger.error("Error c")
                MessageUtils.showMessage("Error encountered reading features: " + e.getMessage());

            }
        } else {
            MessageUtils.showMessage("To use track panning you must first select a single feature track.");
        }


    }

}
