/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
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
