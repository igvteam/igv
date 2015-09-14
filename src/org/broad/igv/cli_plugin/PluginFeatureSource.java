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

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import htsjdk.tribble.Feature;

import javax.xml.bind.annotation.XmlAttribute;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * A feature source which derives its information
 * from a command line cli_plugin
 * User: jacob
 * Date: 2012/05/01
 */
public class PluginFeatureSource<E extends Feature, D extends Feature> extends PluginSource implements FeatureSource{

    private static Logger log = Logger.getLogger(PluginFeatureSource.class);

    @XmlAttribute
    private boolean forbidEmptyOutput;

    @SubtlyImportant
    private PluginFeatureSource(){}

    public PluginFeatureSource(List<String> commands, LinkedHashMap<Argument, Object> arguments, PluginSpecReader.Output outputAttrs, String specPath) {
        this(commands, arguments, outputAttrs, specPath, false);
    }

    public PluginFeatureSource(List<String> commands, LinkedHashMap<Argument, Object> arguments, PluginSpecReader.Output outputAttrs, String specPath, boolean forbidEmptyOutput) {
        super(commands, arguments, outputAttrs, specPath);
        this.forbidEmptyOutput = forbidEmptyOutput;
    }

    @Override
    protected String createTempFile(Track track, Argument argument, String chr, int start, int end, int zoom) throws IOException {
        if(track instanceof AlignmentTrack){
            List<Alignment> alignments = getAlignmentsForRange((AlignmentTrack) track, chr, start, end, zoom);
            return super.createTempFile(alignments, argument);
        }

        FeatureTrack fTrack = (FeatureTrack) track;
        List<Feature> features = fTrack.getFeatures(chr, start, end);
        //Workaround for BEDTools bug, github #88, it can't read an empty file
        if(features.size() == 0 && forbidEmptyOutput){
            features = Arrays.<Feature>asList(new Locus("XXXchr0XXX", 0, 1));
        }
        return super.createTempFile(features, argument);
    }

    /**
     * Perform the actual combination operation between the constituent data
     * sources. This implementation re-runs the operation each call.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws java.io.IOException
     */
    @Override
    public final Iterator<D> getFeatures(String chr, int start, int end) throws IOException {
        return super.getFeatures(chr, start, end, -1);
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    /**
     * If this track has not had it's feature window size set,
     * we use the minimum of the sources
     *
     * @return
     */
    @Override
    public final int getFeatureWindowSize() {
        int featureWindowSize = getMinSizeFromTracks(arguments.values());
        return featureWindowSize;
    }

    private int getMinSizeFromTracks(Iterable tracks) {
        int featureWindowSize = Integer.MAX_VALUE;
        for (Object val : tracks) {
            int tmpSize = Integer.MAX_VALUE;
            if (val instanceof FeatureTrack) {
                FeatureTrack track = ((FeatureTrack) val);
                tmpSize = track.getFeatureWindowSize();
            } else if (val instanceof List) {
                featureWindowSize = getMinSizeFromTracks((List) val);
            }
            featureWindowSize = Math.min(featureWindowSize, tmpSize);
        }
        return featureWindowSize;
    }

    @Override
    public final void setFeatureWindowSize(int size) {
        //no-op
    }


}
