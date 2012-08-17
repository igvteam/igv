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

package org.broad.igv.dev.affective;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

/**
 * Chromosome for the "Affective" project.  A "chromosome" in this project represents part of a calendar day.
 * <p/>
 * Date format is yyyy-MM-dd
 * <p/>
 * e.g. start time
 * Start Time: 2011-04-06 08:57:35 Offset:-04
 *
 * @author Jim Robinson
 * @date 1/21/12
 */
public class AffectiveChromosome implements Chromosome {

    public static DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");

    Date date;
    String name;
    int startTime;  // Start time of this chromosome in seconds.

    int samplingRate = 8;  // # of points per second
    int length;

    public AffectiveChromosome(String dateString) {

        try {
            this.date = dateFormat.parse(dateString);
        } catch (ParseException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        startTime = (AffectiveUtils.START_TIME_HR * 60) * 60;
        length = AffectiveUtils.DAY_LENGTH_HOURS * 60 * 60 * samplingRate;
        name = dateString;
    }


    public int getLength() {
        return length;
    }

    public String getName() {
        return name;
    }

    public List<Cytoband> getCytobands() {
        // Return nothing for now.  This could be used to mark time intervals,  e.g. hours
        return null;
    }

    @Override
    public int getIndex() {
        return 0;
    }

    @Override
    public void setIndex(int ii) {
        //
    }
}
