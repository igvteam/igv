package org.broad.igv.hic.track;

import org.broad.igv.Globals;

/**
 * @author jrobinso
 *         Date: 9/14/12
 *         Time: 8:54 AM
 */
public class HiCFixedGridAxis implements HiCGridAxis {

    int binCount;
    int binSize;
    int igvZoom;

    public HiCFixedGridAxis(int binCount, int binSize) {

        this.binCount = binCount;
        this.binSize = binSize;

        // Compute an approximate igv zoom level
        igvZoom = (int) (Math.log(binCount / 700) / Globals.log2);

    }

    @Override
    public int getGenomicStart(int binNumber) {
        return binNumber * binSize;
    }

    @Override
    public int getGenomicEnd(int binNumber) {
        return binNumber * binSize + binSize;
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
    public int getBinNumberForGenomicPosition(int genomicPosition, int startBin, int endBin) {
        if (endBin < binCount) {
            return (int) (genomicPosition / ((double) binSize));
        } else {
            return -1;
        }
    }

    @Override
    public int getBinNumberForGenomicPosition(int genomicPosition) {
        return (int) (genomicPosition / ((double) binSize));
    }

    @Override
    public int getBinCount() {
        return binCount;
    }
}
