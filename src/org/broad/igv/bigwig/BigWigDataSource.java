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

package org.broad.igv.bigwig;

import org.broad.igv.bbfile.*;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 * @date Jun 18, 2011
 */
public class BigWigDataSource implements DataSource {

    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max);
    WindowFunction windowFunction = WindowFunction.mean;
    Genome genome;
    BBFileReader reader;
    private final BBZoomLevels levels;

    public BigWigDataSource(String path, Genome genome) throws IOException {
        this.genome = genome;

        SeekableStream ss = SeekableStreamFactory.getStreamFor(path);
        reader = new BBFileReader(path, ss);
        levels = reader.getZoomLevels();
    }

    public double getDataMax() {
        return 100;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public double getDataMin() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public List<LocusScore> getSummaryScoresForRange(String chr, int start, int end, int zoom) {
        // Convert zoom level to bp/pixel.
        Chromosome c = genome.getChromosome(chr);
        if (c == null) {
            return null;
        }
        int l = c.getLength();
        double scale = l / (Math.pow(2, zoom) * 700);

        // Lookup closest zoomlevel.  bb zoom level order is opposite of IGVs, proceeds from high to low resolution.
        //
        int bbLevel = -1;
        for (BBZoomLevelHeader zlHeader : levels.getZoomLevelHeaders()) {
            int rl = zlHeader.getReductionLevel();
            if (rl > scale) {
                break;

            }
            else {
                bbLevel = zlHeader.getZoomLevel();
            }
        }

        if (bbLevel >= 0) {
            ArrayList<LocusScore> scores = new ArrayList(1000);
            ZoomLevelIterator zlIter = reader.getZoomLevelIterator(1, chr, start, chr, end, false);
            while (zlIter.hasNext()) {
                ZoomDataRecord rec = zlIter.next();

                // TODO -- check window function
                float mean = (float) (rec.getSumData() / rec.getBasesCovered());

                BasicScore bs = new BasicScore(rec.getChromStart(), rec.getChromEnd(), mean);
                scores.add(bs);
            }
            return scores;

        } else {
            // todo -- compute on the fly
        }


        return null;
    }

    public TrackType getTrackType() {
        return TrackType.OTHER;
    }

    public void setWindowFunction(WindowFunction statType) {
        this.windowFunction = statType;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public void refreshData(long timestamp) {

    }

    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return null;
    }
}
