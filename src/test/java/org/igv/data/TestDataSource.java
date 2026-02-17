package org.igv.data;

import org.igv.feature.LocusScore;
import org.igv.track.TrackType;

import java.util.List;

public class TestDataSource extends AbstractDataSource {


    private int nPts;
    int[] starts;
    int[] ends;
    float[] values;
    String[] probes;

    TestDataSource(int[] starts, int[] ends, float[] values) {
        super(null);
        nPts = starts.length;
        this.starts = starts;
        this.ends = ends;
        this.values = values;
        this.probes = null;
    }

    public TestDataSource() {
        super(null);
        nPts = 10000;
        starts = new int[nPts];
        ends = new int[nPts];
        values = new float[nPts];
        probes = new String[nPts];
        for (int i = 0; i < nPts; i++) {
            starts[i] = i;
            ends[i] = i + 1;
            values[i] = (float) (9.5 + Math.random());
            probes[i] = "probe_" + i;
        }


    }

    protected DataTile getRawData(String chr, int startLocation, int endLocation) {
        return new DataTile(starts, ends, values, probes);
    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public int getLongestFeature(String chr) {
        return 1000;
    }

    public double getMedian(int zoom, String chr) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public double getDataMax() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public double getDataMin() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public TrackType getTrackType() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
