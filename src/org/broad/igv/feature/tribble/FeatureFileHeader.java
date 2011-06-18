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

package org.broad.igv.feature.tribble;

import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: May 22, 2010
 * Time: 7:10:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class FeatureFileHeader {

    /* An object to collection track properties, if specified in the feature file. */
    private TrackProperties trackProperties;

    private TrackType trackType;

    private Set<String> featuresToHide = new HashSet();

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public TrackType getTrackType() {
        return trackType;
    }

    public void setTrackProperties(TrackProperties trackProperties) {
        this.trackProperties = trackProperties;
    }

    public void setTrackType(TrackType trackType) {
        this.trackType = trackType;
    }

    public Set<String> getFeaturesToHide() {
        return featuresToHide;
    }

    public void setFeaturesToHide(Collection<String> featuresToHide) {
        this.featuresToHide.addAll(featuresToHide);
    }
}
