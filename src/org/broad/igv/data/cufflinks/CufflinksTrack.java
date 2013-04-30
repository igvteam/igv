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

package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.Track;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 5:45 PM
 */
public class CufflinksTrack extends DataTrack {

    CufflinksDataSource dataSource;

    public CufflinksTrack(ResourceLocator locator, String id, String name, CufflinksDataSource dataSource) {
        super(locator, id, name);
        this.dataSource = dataSource;
//        setTrackType(TrackType.FPKM);

//        if(isExpDiff(locator.getPath())) {
//            // Make +/- scale symmetic
//            float range = (float) Math.min(5, Math.abs(Math.max(dataSource.getDataMin(), dataSource.getDataMax())));
//            setDataRange(new DataRange(-range, 0, range));
//            setColor(Color.RED);
//            setAltColor(Color.BLUE);
//            setColorScale(new ContinuousColorScale(-range, 0, range, Color.RED, Color.WHITE,  Color.BLUE));
//        }
//
//        else {
//            setDataRange(new DataRange((float) dataSource.getDataMin(), 0, (float) dataSource.getDataMax()));
//        }

    }

    @Override
    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
    }


    // Bit of a hack
    static boolean isExpDiff(String path) {
        return path != null && path.toLowerCase().endsWith("_exp.diff");
    }

    public static void setCufflinksScale(Track inputTrack){
        String path = inputTrack.getResourceLocator() != null ? inputTrack.getResourceLocator().getPath() : null;
        if(isExpDiff(path)) {
            // Make +/- scale symmetic
            float range = Math.min(5, Math.abs(Math.max(inputTrack.getDataRange().getMinimum(), inputTrack.getDataRange().getMaximum())));
            inputTrack.setDataRange(new DataRange(-range, 0, range));
            inputTrack.setColor(Color.RED);
            inputTrack.setAltColor(Color.BLUE);
            inputTrack.setColorScale(new ContinuousColorScale(-range, 0, range, Color.RED, Color.WHITE, Color.BLUE));
        }

        else {
            inputTrack.setDataRange(new DataRange(inputTrack.getDataRange().getMinimum(), 0, inputTrack.getDataRange().getMaximum()));
        }
    }



}
