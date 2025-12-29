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
