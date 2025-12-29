/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.data;

import org.igv.track.TrackProperties;
import org.igv.track.TrackType;

/**
 * A dataset is an in-memory representation of a numerical dataset.  It is used for non-indexed data formats that are
 * read into memory as opposed to read off disk as needed.   The use of objects is eschewed in favor of simple arrays
 * to minimize memory.  
 *
 *
 * @author jrobinso
 */
public interface Dataset {

    public String getName();

    public TrackType getType();

    public TrackProperties getTrackProperties();
    
    public float getDataMin();

    public float getDataMax();

    public String[] getChromosomes();

    public String[] getTrackNames();

    public int[] getStartLocations(String chr);

    public int[] getEndLocations(String chr);

    public String[] getFeatureNames(String chr);

    public float[] getData(String trackName, String chr);

    public boolean isLogNormalized();

    Integer getLongestFeature(String chr);
}
