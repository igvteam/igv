package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.LoadedDataInterval;
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
    }



    @Override
    public LoadedDataInterval<List<LocusScore>>  getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        List<LocusScore> scores =  dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        return new LoadedDataInterval<>(chr, startLocation, endLocation, zoom, scores);
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
