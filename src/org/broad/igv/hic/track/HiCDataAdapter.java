package org.broad.igv.hic.track;

import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;

import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * An adapter class to serve as a bridge between an IGV data source and a HiC track.  HiC tracks differ from
 * IGV tracks in that the coordinate system is based on "bins", each of which can correspond to a variable
 * genomic length.
 *
 * @author jrobinso
 * Date: 9/10/12
 */
public class HiCDataAdapter {

    HiCGridAxis gridAxis;
    DataSource dataSource;

    public HiCDataAdapter(HiCGridAxis gridAxis, DataSource dataSource) {
        this.gridAxis = gridAxis;
        this.dataSource = dataSource;
    }

    public SortedMap<Integer, WeightedSum> getData(String chr, int startBin, int endBin) {

        int zoom = gridAxis.getIGVZoom();
        int gStart = gridAxis.getGenomicStart(startBin);
        int gEnd = gridAxis.getGenomicEnd(endBin);

        List<LocusScore> scores = dataSource.getSummaryScoresForRange(chr, gStart, gEnd, zoom);

        // Compute weighted average of signal over each bin
        SortedMap<Integer, WeightedSum> dataBinMap =  new TreeMap<Integer, WeightedSum>();

        for (LocusScore locusScore : scores) {

            int bs = gridAxis.getBinNumberForGenomicPosition(locusScore.getStart());
            int be = gridAxis.getBinNumberForGenomicPosition(locusScore.getEnd());

            if(bs > endBin) {
                break;
            }
            else if(be < startBin) {
                continue;
            }

            for (int b = bs; b <= be; b++) {

                int bStart = gridAxis.getGenomicStart(b);
                int bEnd = gridAxis.getGenomicEnd(b);
                WeightedSum dataBin = dataBinMap.get(b);
                if(dataBin == null) {
                    dataBin = new WeightedSum(b, bStart, bEnd);
                    dataBinMap.put(b, dataBin);
                }
                dataBin.addScore(locusScore);

            }
        }

        return dataBinMap;

    }

    public static class WeightedSum {
        int binNumber;
        int nPts = 0;
        double weightedSum = 0;
        int genomicStart;
        int genomicEnd;

        private WeightedSum(int binNumber, int genomicStart, int genomicEnd) {
            this.binNumber = binNumber;
            this.genomicStart = genomicStart;
            this.genomicEnd = genomicEnd;
        }

        public int getBinNumber() {
            return binNumber;
        }

        void addScore(LocusScore ls) {
            if(ls.getStart() >= genomicEnd || ls.getEnd() < genomicStart) return;

            double weight = ((double) (Math.min(genomicEnd, ls.getEnd()) - Math.max(genomicStart, ls.getStart()))) /
                    (genomicEnd - genomicStart);
            weightedSum += weight * ls.getScore();
            nPts++;
        }


        public float weightedAverage() {
            return nPts == 0 ? 0 : (float) (weightedSum / nPts);
        }


    }


}
