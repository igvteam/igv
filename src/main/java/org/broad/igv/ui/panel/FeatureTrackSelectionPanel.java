package org.broad.igv.ui.panel;

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.util.Collection;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/14/12
 * Time: 8:34 AM
 * To change this template use File | Settings | File Templates.
 */
public class FeatureTrackSelectionPanel extends JPanel {


    public FeatureTrackSelectionPanel() {
        init();
    }

    void init() {

        if(IGV.hasInstance()) {

            LayoutManager lm = new BoxLayout(this, BoxLayout.Y_AXIS);
            this.setLayout(lm);

            Collection<Track> tracks = IGV.getInstance().getAllTracks();
            for(Track t : tracks) {
                if(t instanceof FeatureTrack) {

                    JCheckBox cb = new JCheckBox(t.getName());
                    this.add(cb);


                }
            }
        }

    }

}
