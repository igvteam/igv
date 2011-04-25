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

package org.broad.igv.peaks;

import org.apache.log4j.Logger;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;


import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;

/**
 * @author jrobinso
 * @date Apr 23, 2011
 */
public class PeakTrackMenu extends JPopupMenu {

    private static Logger log = Logger.getLogger(PeakTrackMenu.class);
    private PeakTrack track;

    public PeakTrackMenu(PeakTrack t) {
        this.track = t;

        //Title
        JLabel popupTitle = new JLabel("<html><b>" + track.getName(), JLabel.LEFT);
        Font newFont = getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        add(popupTitle);

        //Change Track Settings
        addDisplayModeItems();
        addColorByItems();
        addRendererItems();

        addSeparator();
        JMenuItem item = new JMenuItem("Open control panel");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                PeakTrack.openControlDialog();
            }
        });
        item.setEnabled(!PeakTrack.controlDialogIsOpen());
        add(item);

        addSeparator();
        add(TrackMenuUtils.getRemoveMenuItem(Arrays.asList(new Track[]{track})));

    }

    public void addDisplayModeItems() {
        addSeparator();
         ButtonGroup group = new ButtonGroup();

        Track.DisplayMode displayMode = track.getDisplayMode();

        JRadioButtonMenuItem m1 = new JRadioButtonMenuItem("Collapse");
        m1.setSelected(displayMode == Track.DisplayMode.COLLAPSED);
        m1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.COLLAPSED);
                IGV.getInstance().repaint();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("Expand");
        m3.setSelected(displayMode == Track.DisplayMode.EXPANDED);
        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.EXPANDED);
                IGV.getInstance().repaint();
            }
        });

        add(m1);
        add(m3);
        group.add(m1);
        group.add(m3);

    }

    public void addColorByItems() {


        addSeparator();
        add(new JLabel("Color by"));

        ButtonGroup group = new ButtonGroup();

        JRadioButtonMenuItem m1 = new JRadioButtonMenuItem("Score");
        m1.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.SCORE);
        m1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                PeakTrack.setShadeOption(PeakTrack.ColorOption.SCORE);
                IGV.getInstance().repaint();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("Fold change");
        m3.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.FOLD_CHANGE);
        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                PeakTrack.setShadeOption(PeakTrack.ColorOption.FOLD_CHANGE);
                IGV.getInstance().repaint();
            }
        });

        add(m1);
        add(m3);
        group.add(m1);
        group.add(m3);

    }

        public void addRendererItems() {


        addSeparator();
        add(new JLabel("Display as"));

        ButtonGroup group = new ButtonGroup();

        JRadioButtonMenuItem m1 = new JRadioButtonMenuItem("Features");
        m1.setSelected(track.getRenderOption() == PeakTrack.RenderOption.FEATURE);
        m1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                PeakTrack.setRenderOption(PeakTrack.RenderOption.FEATURE);
                IGV.getInstance().repaint();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("Barchart");
        m3.setSelected(track.getRenderOption() == PeakTrack.RenderOption.CHART);
        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                PeakTrack.setRenderOption(PeakTrack.RenderOption.CHART);
                IGV.getInstance().repaint();
            }
        });

        add(m1);
        add(m3);
        group.add(m1);
        group.add(m3);

    }
}
