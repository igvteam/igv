/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.maf.conservation;

import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.DataSource;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.LRUCache;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public class OmegaDataSource extends AbstractDataSource implements DataSource {

    //todo this is where we should put this data, but it's not there
    File rootDir = new File("/resources/omega/12mer");

    public OmegaDataSource(Genome genome) {
        super(genome);
    }

    public double getDataMax() {
        return 1;
    }

    public double getDataMin() {
        return 0;
    }

    public int getLongestFeature(String chr) {
        return 1;
    }


    public TrackType getTrackType() {
        return TrackType.OTHER;
    }

    public void setWindowFunction(WindowFunction statType) {
        // ignored for now
    }

    public boolean isLogNormalized() {
        return false;
    }

    public void refreshData(long timestamp) {
        // ignore
    }

    public WindowFunction getWindowFunction() {
        return null;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return new ArrayList();
    }

    public String getFileNameForPosition(String chr, int pos) {
        //chr1_2000000-2999999.omega
        int mb = (pos / 1000000) * 1000000;
        int end = mb + 999999;
        return chr + "_" + mb + "-" + end + ".omega";
    }


    LRUCache<String, DataTile> tileCache = new LRUCache(this, 3);

    @Override
    public DataTile getRawData(String chr, int startLocation, int endLocation) {

        int bStart = startLocation / 1000000;
        int bEnd = endLocation / 1000000;

        String key = chr + bStart + "_" + bEnd;
        if (tileCache.containsKey(key)) {
            return tileCache.get(key);
        }

        File chrDir = new File(rootDir, chr.substring(3));
        List<File> files = new ArrayList<File>(bEnd - bStart + 1);
        for (int i = bStart; i <= bEnd; i++) {

            int start = i * 1000000;
            int end = start + 999999;
            File f = new File(chrDir, chr + "_" + start + "-" + end + ".omega");
            if (f.exists()) {
                files.add(f);
            }

        }

        if (files.isEmpty()) {
            return null;
        }

        int sz = endLocation - startLocation;
        int[] starts = new int[sz];
        float[] data = new float[sz];

        OmegaFileParser parser = new OmegaFileParser();

        int offset = 0;
        for (File f : files) {
            parser.getDataAsArrays(f, starts, data, offset, startLocation, endLocation);
            offset += 1000000;
        }
        DataTile dt = new DataTile(starts, null, data, null);
        tileCache.put(key, dt);
        return dt;

    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }

}
