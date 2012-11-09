/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.hic.track;

import org.broad.igv.Globals;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * @author jrobinso
 *         Date: 9/14/12
 *         Time: 8:49 AM
 */
public class HiCFragmentAxis implements HiCGridAxis {

    double averageBinSize;
    int igvZoom;
    int[] sites;
    int chrLength;


    /**
     * @param sites  ordered by start position.  Its assumed bins are contiguous, no gaps and no overlap.
     * @param length
     */
    public HiCFragmentAxis(int[] sites, int length) {

        this.sites = sites;

        this.chrLength = length;

        averageBinSize = chrLength / (sites.length + 1);

        // Compute an approximate igv zoom level

        igvZoom = (int) (Math.log((chrLength / 700) / averageBinSize) / Globals.log2);


    }


    @Override
    public int getGenomicStart(int binNumber) {

        return binNumber == 0 ? 0 : sites[binNumber - 1];

    }

    @Override
    public int getGenomicEnd(int binNumber) {
        return binNumber < sites.length ? sites[binNumber] : chrLength;
    }

    @Override
    public int getGenomicMid(int binNumber) {
        return (getGenomicStart(binNumber) + getGenomicEnd(binNumber)) / 2;
    }


    @Override
    public int getIGVZoom() {
        return igvZoom;
    }


    /**
     * Return fragment (bin) that this position lies on.  Fragment 0 means position < sites[0].
     * Fragment 1 means position >= sites[0] and < sites[1].
     *
     * @param position The genome position to search for within that array
     * @return The fragment location such that position >= sites[retVal-1] and position <  sites[retVal]
     */
    @Override
    public int getBinNumberForGenomicPosition(int position) {
        int lo = 0;
        int hi = sites.length - 1;
        while (lo <= hi) {
            // Base case - found range
            int mid = lo + (hi - lo) / 2;

            if (position > sites[mid]) lo = mid + 1;
            else if (position < sites[mid]) hi = mid - 1;
            else return mid + 1;
        }
        return lo;
    }

    @Override
    public int getBinNumberForFragment(int fragment) {
        if (fragment <= sites.length) {
            return fragment;
        } else {
            throw new RuntimeException("Fragment: " + fragment + " is out of range");
        }
    }

    @Override
    public int getBinCount() {
        return sites.length + 1;
    }


}
