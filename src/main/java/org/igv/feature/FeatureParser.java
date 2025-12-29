/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.feature;

//~--- JDK imports ------------------------------------------------------------


import org.igv.feature.genome.Genome;
import org.igv.track.FeatureTrack;
import org.igv.track.TrackProperties;
import org.igv.util.ResourceLocator;

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
