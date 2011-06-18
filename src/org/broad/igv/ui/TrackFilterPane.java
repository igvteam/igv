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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui;

import org.broad.igv.track.Track;
import org.broad.igv.ui.util.FilterComponent;
import org.broad.igv.ui.util.FilterPane;
import org.broad.igv.util.Filter;
import org.broad.igv.util.FilterElement;

import java.awt.*;
import java.util.List;

/**
 * @author eflakes
 */
public class TrackFilterPane extends FilterPane {

    private boolean matchAll = true;

    public TrackFilterPane(List<String> items, String itemListLabel, TrackFilter filter) {
        super(items, itemListLabel, filter);
    }

    public Filter createNewFilter() {
        String name = "";
        return new TrackFilter(name, null);
    }

    @Override
    public TrackFilter getFilter() {
        return (TrackFilter) super.getFilter();
    }

    public void clearTracks() {
        getFilter().clearTracks();
    }

    public void addTracks(List<Track> tracks) {
        getFilter().addTracks(tracks);
    }

    public FilterComponent createFilterComponent(FilterPane filterPane,
                                                 String itemListLabel,
                                                 List<String> itemList,
                                                 FilterElement element) {

        return new TrackFilterComponent((TrackFilterPane) filterPane, itemListLabel, itemList, element);
    }

    public void setMatchAll(boolean value) {
        matchAll = value;
    }

    public boolean getMatchAll() {
        return matchAll;
    }

    public void applyFilterMatching() {
        Component[] filterComponents = getComponents();
        for (Component filterComponent : filterComponents) {
            ((TrackFilterComponent) filterComponent).setMatchAll(matchAll);
        }
    }

    /**
     * Only for use with more() method.
     */
    public void addDefaultFilterComponent() {

        if (itemListLabel == null) {
            itemListLabel = "When";
        }
        add(new TrackFilterComponent(this, itemListLabel, itemList, null));
    }

}
