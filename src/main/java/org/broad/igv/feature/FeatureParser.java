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

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome);


    TrackProperties getTrackProperties();

}
