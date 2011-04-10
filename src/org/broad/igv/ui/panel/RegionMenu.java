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

package org.broad.igv.ui.panel;

import com.jidesoft.swing.JidePopupMenu;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Set;

/**
 * @author jrobinso
 * @date Sep 18, 2010
 */
public class RegionMenu extends JidePopupMenu {

    RegionOfInterest roi;
    ReferenceFrame frame;

    public RegionMenu(final RegionOfInterest roi, final ReferenceFrame frame) {
        this(roi, frame, null);
    }

    public RegionMenu(final RegionOfInterest roi, final ReferenceFrame frame, String title) {

        this.roi = roi;
        this.frame = frame;

        if (title != null) {
            JLabel heading = new JLabel("<html>&nbsp;<b>" + title);
            //heading.setFont(UIConstants.boldFont);
            add(heading);
            addSeparator();
        }


        Set<TrackType> loadedTypes = IGV.getInstance().getTrackManager().getLoadedTypes();

        if (loadedTypes.contains(TrackType.COPY_NUMBER) ||
                loadedTypes.contains(TrackType.ALLELE_SPECIFIC_COPY_NUMBER) ||
                loadedTypes.contains(TrackType.CNV)) {
            JMenuItem item = new JMenuItem("Sort by amplification");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().getTrackManager().sortByRegionScore(roi, RegionScoreType.AMPLIFICATION, frame);

                }
            });
            add(item);


            item = new JMenuItem("Sort by deletion");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().getTrackManager().sortByRegionScore(roi, RegionScoreType.DELETION, frame);
                    IGV.getInstance().getContentPane().repaint();
                }
            });
            add(item);

            item = new JMenuItem("Sort by breakpoint amplitudes");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().getTrackManager().sortByRegionScore(roi, RegionScoreType.FLUX, frame);
                    IGV.getInstance().getContentPane().repaint();
                }
            });
            add(item);
        }

        if (loadedTypes.contains(TrackType.GENE_EXPRESSION)) {
            JMenuItem item = new JMenuItem("Sort by expression");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().getTrackManager().sortByRegionScore(roi, RegionScoreType.EXPRESSION, frame);
                    IGV.getInstance().getContentPane().repaint();

                }
            });

            add(item);
        }

        if (loadedTypes.contains(TrackType.MUTATION)) {
            JMenuItem item = new JMenuItem("Sort by muation count");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().getTrackManager().sortByRegionScore(roi, RegionScoreType.MUTATION_COUNT, frame);
                    IGV.getInstance().getContentPane().repaint();

                }
            });

            add(item);
        }

        JMenuItem item = new JMenuItem("Sort by value");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                IGV.getInstance().getTrackManager().sortByRegionScore(roi, RegionScoreType.SCORE, frame);
                IGV.getInstance().getContentPane().repaint();

            }
        });
        add(item);

    }


}
