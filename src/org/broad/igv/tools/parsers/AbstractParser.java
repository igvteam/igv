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

package org.broad.igv.tools.parsers;

import org.broad.igv.track.TrackType;

import java.util.Hashtable;

/**
 * @author jrobinso
 */
public abstract class AbstractParser {

    private static Hashtable<String, String> chromosomeNames = new Hashtable();

    private DataConsumer dataConsumer;
    String trackLine = null;
    private TrackType trackType = TrackType.COPY_NUMBER;
    private String[] headings;

    public AbstractParser(DataConsumer dataConsumer) {
        this.dataConsumer = dataConsumer;
    }

    /**
     * @return the dataConsumer
     */
    public DataConsumer getDataConsumer() {
        return dataConsumer;
    }

    protected void setTrackParameters() {
        dataConsumer.setTrackParameters(trackType, trackLine, headings);
    }

    /**
     * Note:  This is an exact copy of the method in ExpressionFileParser.  Refactor to merge these
     * two parsers, or share a common base class.
     *
     * @param comment
     */
    void parseComment(String comment) {
        String tmp = comment.substring(1, comment.length());
        if (tmp.startsWith("track")) {
            trackLine = tmp;
            System.out.println("track line = " + trackLine);
        } else {
            String[] tokens = tmp.split("=");
            if (tokens.length != 2) {
                return;
            }
            String key = tokens[0].trim().toLowerCase();
            if (key.equals("type")) {
                try {
                    setTrackType(TrackType.valueOf(tokens[1].trim().toUpperCase()));
                } catch (Exception exception) {
                }
            }
        }
    }

    /**
     * @param headings the headings to set
     */
    public void setHeadings(String[] headings) {
        this.headings = headings;
    }

    /**
     * @return the headings
     */
    public String[] getHeadings() {
        return headings;
    }

    /**
     * @return the trackType
     */
    public TrackType getTrackType() {
        return trackType;
    }

    /**
     * @param trackType the trackType to set
     */
    public void setTrackType(TrackType trackType) {
        this.trackType = trackType;
    }

}
