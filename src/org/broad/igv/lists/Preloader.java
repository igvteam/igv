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

package org.broad.igv.lists;

import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.Collection;

/**
 * This class added to preload data when using gene lists.   Its actually more general than that,  but its motivation
 * stems from the need to provide a wait cursor when loading gene lists.
 *
 * @author jrobinso
 * @date Mar 22, 2011
 */
public class Preloader {

    public static synchronized void preload() {

        Collection<Track> trackList = IGV.getInstance().getAllTracks();
        int flankingRegion = 1; //PreferenceManager.getInstance().getAsInt(PreferenceManager.FLANKING_REGION) + 1;
        String genomeId = GenomeManager.getInstance().getGenomeId();
        for (ReferenceFrame frame : FrameManager.getFrames()) {
            Locus locus = frame.getInitialLocus();
            if (locus != null) {
                for (Track track : trackList) {
                    if (track == null) continue;
                    if (track.isVisible()) {
                        if (track instanceof DataTrack) {
                            DataTrack dt = (DataTrack) track;
                            RenderContext context = new RenderContextImpl(null, null, frame, null);
                            int start = Math.max(0, locus.getStart() - flankingRegion);
                            int end = locus.getEnd() + flankingRegion;
                            dt.load(context, locus.getChr(), start, end, frame.getZoom());
                        }
                    }
                }
            }
        }
    }
}

