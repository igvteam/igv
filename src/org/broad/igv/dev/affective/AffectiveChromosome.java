package org.broad.igv.dev.affective;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.List;

/**
 * Chromosome for the "Affective" project.  A "chromosome" in this project represents part of a calendar day.
 *
 * e.g. start time
 *   Start Time: 2011-04-06 08:57:35 Offset:-04
 * @author Jim Robinson
 * @date 1/21/12
 */
public class AffectiveChromosome implements Chromosome {

    Date date;
    String name;
    //int startTime;
    int pointSpacing = 8;  // # of points per second
    int length;

    public AffectiveChromosome(Date date) {

        this.date = date;
        //startTime = AffectiveUtils.START_TIME;
        length = AffectiveUtils.DAY_LENGTH_HOURS * 60 * 60 * pointSpacing;

        name = (new SimpleDateFormat("EEE, MMM d")).format(date);
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
}
