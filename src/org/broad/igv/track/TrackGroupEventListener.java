package org.broad.igv.track;

import java.util.EventListener;

/**
 * @author Jim Robinson
 * @date 10/28/11
 */
public interface TrackGroupEventListener extends EventListener {

    public  void onTrackGroupEvent(TrackGroupEvent e);
}
