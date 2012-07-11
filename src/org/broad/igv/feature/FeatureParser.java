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


/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/

package org.broad.igv.feature;

//~--- JDK imports ------------------------------------------------------------


import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;

/**
 * @author jrobinso
 */
public interface FeatureParser {


    List<FeatureTrack> loadTracks(ResourceLocator locator, Genome genome) throws IOException;

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    List<org.broad.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome);


    public TrackProperties getTrackProperties();

}
