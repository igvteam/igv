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
    int[] sites;

    public HiCFixedGridAxis(int binCount, int binSize, int [] sites) {

        this.binCount = binCount;
        this.binSize = binSize;
        this.sites = sites;

        // Compute an approximate igv zoom level
        igvZoom = Math.max(0, (int) (Math.log(binCount / 700) / Globals.log2));

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
    public int getGenomicMid(int binNumber) {
        return binNumber * binSize + binSize / 2;
    }

    @Override
    public int getIGVZoom() {
        return igvZoom;
    }

    @Override
    public int getBinNumberForGenomicPosition(int genomicPosition) {
        return (int) (genomicPosition / ((double) binSize));
    }

    @Override
    public int getBinNumberForFragment(int fragment) {

        if (fragment < sites.length) {
            int genomicPosition = sites[fragment];
            return getBinNumberForGenomicPosition(genomicPosition);
        }
        throw new RuntimeException("Fragment: " + fragment + " is out of range");
    }

    @Override
    public int getBinCount() {
        return binCount;
    }

}
