package org.broad.igv.ui.event;

import java.util.EventListener;

/**
 * @author Jim Robinson
 * @date 12/2/11
 */
public interface AlignmentTrackEventListener extends EventListener {

    public  void onAlignmentTrackEvent(AlignmentTrackEvent e);

}
