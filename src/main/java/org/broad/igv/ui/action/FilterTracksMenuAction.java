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
package org.broad.igv.ui.action;

import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.TrackFilterEvent;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilterDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.TrackFilter;

import java.awt.event.ActionEvent;
import java.util.List;

/**
 * @author jrobinso
 */
public class FilterTracksMenuAction extends MenuAction {

    //static Logger log = LogManager.getLogger(FilterTracksMenuAction.class);
    IGV igv;


    public FilterTracksMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        List<String> uniqueAttributeKeys = AttributeManager.getInstance().getAttributeNames();

        // Sort the attribute keys if we have any
        if (uniqueAttributeKeys != null) {
            //Collections.sort(uniqueAttributeKeys, AttributeManager.getInstance().getAttributeComparator());
        } else // If we have no attribute we can't display the
            // track filter dialog so say so and return
            if (uniqueAttributeKeys == null || uniqueAttributeKeys.isEmpty()) {

                MessageUtils.showMessage("No attributes found to use in a filter");
                return;
            }

        TrackFilter trackFilter = IGV.getInstance().getSession().getFilter();
        TrackFilterDialog dialog = new TrackFilterDialog(igv.getMainFrame(), "Filter Tracks", trackFilter);
        dialog.setVisible(true);

        if (!dialog.isCancelled()) {
            trackFilter = dialog.getFilter();
            filterTracks(trackFilter);
            igv.repaint();
        }

    }



    private void filterTracks(TrackFilter trackFilter) {

        if (trackFilter == null) {
            List<Track> tracks = IGV.getInstance().getAllTracks();
            IGV.getInstance().getSession().setFilter(null);
            for (Track track : tracks) {
                track.setVisible(true);
            }
        } else {
            IGV.getInstance().getSession().setFilter(trackFilter);
            trackFilter.evaluate();
        }

        for(Track track : IGV.getInstance().getAllTracks()) {
            track.filterSamples(trackFilter);
        }
    }

}
