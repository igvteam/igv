package org.broad.igv.hic.track;

/**
 * @author jrobinso
 *         Date: 9/14/12
 *         Time: 8:54 AM
 */
public interface HiCGridAxis {

    int getGenomicStart(int binNumber);

    int getGenomicEnd(int binNumber);

    int getGenomicMid(int binNumber);

    int getIGVZoom();

    int getBinCount();

    int getBinNumberForGenomicPosition(int genomePosition);

    int getBinNumberForFragment(int fragmentX);
}
