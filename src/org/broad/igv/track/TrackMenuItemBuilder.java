/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

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
