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

package org.broad.igv.data.seg;

import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ResourceLocator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class SegmentedBinaryDataSet implements SegmentedDataSet {

    SegmentedBinaryReader reader;

    private long lastRefreshTime = 0;

    /**
     * Assume data is non-log value until suggested otherwise by the precense
     * of negative numbers.  TODO This is a fragile assumption, the user should
     * input this information directly.
     */
    Boolean logNormalized = null;

    private List<String> sampleNames = null;

    /**
     * Map of [heading ->  [chr -> [list of segmentsCache]]]
     */
    Map<String, Map<String, List<LocusScore>>> segmentsCache = new HashMap();

    Map<String, SegmentedChromosomeData> chrData = new HashMap();

    TrackType type = TrackType.COPY_NUMBER;

    public SegmentedBinaryDataSet(ResourceLocator locator) {

        reader = new SegmentedBinaryLocalReader(locator.getPath());

        try {
            type = TrackType.valueOf(reader.getStringAttribute("type"));
        }
        catch (Exception e) {
            //ignore;
        }

        try {
            String logNormalizedString = reader.getStringAttribute("logNormalized");
            if (logNormalizedString != null) {
                logNormalized = logNormalizedString.equals("true");
            }

        } catch (Exception exception) {
            // ignore
        }


    }


    // TODO -- the synchronized keyword can be removed with the async whole genome load

    public List<LocusScore> getSegments(String heading, String chr) {

        Map<String, List<LocusScore>> chrSegments = segmentsCache.get(heading);
        if (chrSegments == null) {
            chrSegments = new HashMap();
            segmentsCache.put(heading, chrSegments);
        }

        List<LocusScore> segments = chrSegments.get(chr);

        if (segments == null) {

            SegmentedChromosomeData cd = chrData.get(chr);
            if (cd == null) {
                cd = reader.getChromosomeData(chr);
                chrData.put(chr, cd);
            }

            int[] startLocations = cd.getStartLocations(heading);
            int[] endLocations = cd.getEndLocations(heading);
            float[] values = cd.getValues(heading);

            if (startLocations == null || startLocations.length == 0) {
                return null;
            }

            if (logNormalized == null) {
                logNormalized = false;
                for (int i = 0; i <
                        values.length; i++) {
                    if (values[i] < 0) {
                        logNormalized = true;
                        break;
                    }

                }
            }
            assert (startLocations.length == endLocations.length);
            assert (endLocations.length == values.length);

            segments = new ArrayList(startLocations.length);
            for (int i = 0; i <
                    startLocations.length; i++) {
                segments.add(new Segment(startLocations[i], endLocations[i], values[i]));
                chrSegments.put(chr, segments);
            }

        }
        return segments;
    }

    public List<String> getSampleNames() {
        if (sampleNames == null) {
            sampleNames = reader.getSampleNames();
        }

        return sampleNames;
    }

    public TrackType getType() {
        return type;
    }

    // TODO -- check for negative numbers

    public boolean isLogNormalized() {
        return logNormalized == null ? true : logNormalized;
    }

    public double getDataMax(String chr) {
        return 3;
    }

    public double getDataMin(String chr) {
        return -3;
    }

    public synchronized List<LocusScore> getWholeGenomeScores(String heading) {
        return getSegments(heading, Globals.CHR_ALL);

    }

    public static class ChromosomeChunk {
    }

}
