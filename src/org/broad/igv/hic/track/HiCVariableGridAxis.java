package org.broad.igv.hic.track;

import com.sun.tools.internal.xjc.reader.xmlschema.BindGreen;
import org.broad.igv.Globals;
import sun.tools.tree.BinaryShiftExpression;

/**
 * @author jrobinso
 *         Date: 9/14/12
 *         Time: 8:49 AM
 */
public class HiCVariableGridAxis implements HiCGridAxis {

    Bin[] bins;
    double averageBinSize;
    int igvZoom;

    /**
     *
     * @param bins ordered by start position.  Its assumed bins are contiguous, no gaps and no overlap.
     */
    public HiCVariableGridAxis(Bin[] bins) {

        this.bins = bins;

        Bin lastBin = bins[bins.length - 1];
        double chrLength = lastBin.start + lastBin.width;

        averageBinSize = bins.length == 0 ? 0 : (chrLength / bins.length);

        // Compute an approximate igv zoom level

        igvZoom =  (int) (Math.log((chrLength / 700) / averageBinSize) / Globals.log2);


    }

    @Override
    public int getGenomicStart(int binNumber) {

        return bins[binNumber].start;

    }

    @Override
    public int getGenomicEnd(int binNumber) {
        Bin b = bins[binNumber];
        return b.start + b.width;
    }

    /**
     * Return the resolution in base-pairs / pixel.  This is a representative value.
     *
     * @return
     */
    @Override
    public double getResolution() {
        return averageBinSize;
    }

    @Override
    public int getIGVZoom() {
        return igvZoom;
    }

    /**
     * Return the bin number containing the genomic position.  The bin is restrained by startBin and endBin.  If
     * no bins in this range contain the position return -1  (this should not happen).
     */
    @Override
    public int getBinNumberForGenomicPosition(int start, int startBin, int endBin) {
        for(int b = startBin; b <= endBin; b++) {
            Bin bin = bins[b];
            if(start >= bin.start && start < (bin.start + bin.width)) {
                return b;
            }
        }
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public static class Bin {

        public Bin(int start, int width) {
            this.start = start;
            this.width = width;
        }

        int start;
        int width;
    }

}
