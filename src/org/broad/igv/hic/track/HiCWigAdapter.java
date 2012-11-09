package org.broad.igv.hic.track;

import org.broad.igv.data.BasicScore;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.data.WiggleParser;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.hic.HiC;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 *         Date: 11/8/12
 *         Time: 10:30 AM
 */
public class HiCWigAdapter extends HiCDataAdapter {

    WiggleDataset dataset;
    String trackName;
    Map<String, List<LocusScore>> locusScoreMap = new HashMap<String, List<LocusScore>>();
    private Color color;
    private Color altColor;
    DataRange dataRange;

    public HiCWigAdapter(HiC hic, String path) {
        super(hic);
        init(path);
    }


    private void init(String path) {

        dataset = (new WiggleParser(new ResourceLocator(path), null)).parse();

        trackName = dataset.getTrackNames()[0];

        TrackProperties properties = dataset.getTrackProperties();

        color = properties.getColor();
        if (color == null) color = Color.blue.darker();
        altColor = properties.getAltColor();

        float min = properties.getMinValue();
        float max = properties.getMaxValue();
        float mid = properties.getMidValue();
        if (Float.isNaN(min) || Float.isNaN(max)) {
            mid = 0;
            min = dataset.getPercent10();
            max = dataset.getPercent90();
            if(min > 0 && max > 0) min = 0;
            else if(min < 0 && max < 0) max = 0;


        } else {
            if (Float.isNaN(mid)) {
                if (min >= 0) {
                    mid = Math.max(min, 0);
                } else {
                    mid = Math.min(max, 0);
                }
            }
        }

        dataRange = new DataRange(min, mid, max);
        if (properties.isLogScale()) {
            dataRange.setType(DataRange.Type.LOG);
        }

    }


    protected java.util.List<LocusScore> getLocusScores(String chr, int zoom, int gStart, int gEnd) {

        if(!chr.startsWith("chr")) chr = "chr" + chr;

        List<LocusScore> scores = locusScoreMap.get(chr);
        if (scores == null) {
            int[] startLocations = dataset.getStartLocations(chr);
            int[] endLocations = dataset.getEndLocations(chr);
            float[] values = dataset.getData(trackName, chr);
            scores = new ArrayList<LocusScore>(values.length);
            for (int i = 0; i < values.length; i++) {
                BasicScore bs = new BasicScore(startLocations[i], endLocations[i], values[i]);
                scores.add(bs);
            }
            locusScoreMap.put(chr, scores);

        }
        return scores;
    }

    @Override
    public String getName() {
        return trackName;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Color getColor() {
        return color;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean isLogScale() {
        return dataRange == null ? false : dataRange.isLog();  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Color getAltColor() {
        return altColor;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public DataRange getDataRange() {
        return dataRange;  //To change body of implemented methods use File | Settings | File Templates.
    }

}
