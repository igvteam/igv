/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tools.parsers;

import org.broad.igv.track.TrackType;

/**
 * @author jrobinso
 */
public interface DataConsumer {

    public void setType(String type);

    public void addData(String chr, int start, int end, float[] data, String name);

    public void parsingComplete();

    public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames);

    void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames, boolean b);

    /**
     * Set a tolerance for "sortedness" of the data.  A start position can be less than
     * the immediately previous start position by this amount.  This is needed for
     * chip-seq processing where the start position of an alignment can be artificially
     * extended post sorting.
     */
    public void setSortTolerance(int tolerance);

    public void setAttribute(String key, String value);

}
