/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.DataTrack;
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
    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
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
