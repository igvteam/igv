package org.igv.track;

import javax.swing.*;
import java.util.Collection;

/**
 * Interface to implement when adding a context menu item for track context menus
 *
 * @author jacob
 * @date 2012-Dec-21
 * @api
 */
public interface TrackMenuItemBuilder {

    /**
     * @param selectedTracks The tracks currently selected when the user clicks
     * @param te
     * @return A {@code JMenuItem} to be displayed in the menu,
     *         or {@code null} if not applicable to the selectedTracks
     */
    JMenuItem build(Collection<Track> selectedTracks, TrackClickEvent te);
}
