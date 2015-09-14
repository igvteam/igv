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

package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.feature.tribble.MUTCodec;
import org.broad.igv.track.*;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Parses a mutation file, such as ".mut" or ".maf" (mutation annotation file)
 *
 * @author jrobinso
 */
public class MutationTrackLoader {

    private static Logger log = Logger.getLogger(MutationTrackLoader.class);
    private ResourceLocator locator = null;
    private Genome genome = null;
    MUTCodec codec;

    public static boolean isMutationAnnotationFile(ResourceLocator locator) throws IOException {
        return MUTCodec.isMutationAnnotationFile(locator);
    }

    public List<FeatureTrack> loadMutationTracks(ResourceLocator locator, Genome genome) throws IOException, TribbleIndexNotFoundException {

        this.locator = locator;
        this.genome = genome;

        boolean indexed = isIndexed(locator, genome);

        List<FeatureTrack> tracks = new ArrayList<FeatureTrack>();

        if (indexed) {
            String[] samples = getCodec().getSamples();
            MutationFeatureSource.MutationDataManager dataManager = new MutationFeatureSource.MutationDataManager(locator, genome);
            for (String sampleId : samples) {
                String id = locator.getPath() + "_" + sampleId;
                FeatureSource<Mutation> featureSource = new MutationFeatureSource(sampleId, dataManager);
                MutationTrack track = new MutationTrack(locator, id, featureSource);
                tracks.add(track);
                track.setName(sampleId);
            }

        } else {
            Map<String, List<htsjdk.tribble.Feature>> features = loadMutations();
            for (String sampleId : features.keySet()) {
                String id = locator.getPath() + "_" + sampleId;
                MutationTrack track = new MutationTrack(locator, id, new FeatureCollectionSource(features.get(sampleId), genome));
                tracks.add(track);
                track.setName(sampleId);
            }
        }

        for (FeatureTrack track : tracks) {
            track.setSquishedRowHeight(5);
            track.setExpandedRowHeight(15);
            track.setHeight(15);
            // Override default minimum height (10 for feature tracks).
            track.setMinimumHeight(0);
        }
        //Just to make sure we have no memory
        this.locator = null;
        this.genome = null;
        return tracks;
    }


    private MUTCodec getCodec() {
        if (codec == null) codec = new MUTCodec(locator.getPath(), genome);
        return codec;
    }

    /**
     * Test to see if a usable index exists.  In addition to the index, mutation files have an additional requirement
     * that samples be specified in a header directive.
     *
     * @param locator
     * @return
     */
    private boolean isIndexed(ResourceLocator locator, Genome genome) {
        if (!TrackLoader.isIndexed(locator, genome)) return false;

        try {
            String[] samples = getCodec().getSamples();
            return samples != null && samples.length > 0;
        } catch (Exception e) {
            log.error("Error creating codec for: " + locator.getPath(), e);
            return false;
        }

    }

    /**
     * Return a map of runId -> list of mutation objects.   The "runId" field
     * is the track identifier (name) for mutation files.
     *
     * @return
     */
    private Map<String, List<htsjdk.tribble.Feature>> loadMutations() {
        AsciiLineReader reader = null;
        String nextLine = null;

        try {

            if (codec == null) codec = new MUTCodec(locator.getPath(), genome);

            Map<String, List<htsjdk.tribble.Feature>> mutationMap = new LinkedHashMap();

            reader = ParsingUtils.openAsciiReader(locator);

            // Skip header - handled in codec
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;
                else break;

            }

            while ((nextLine = reader.readLine()) != null) {

                Mutation mut = codec.decode(nextLine);

                if (mut != null) {
                    String sampleId = mut.getSampleId();
                    List<htsjdk.tribble.Feature> features = mutationMap.get(sampleId);
                    if (features == null) {
                        features = new ArrayList<Feature>();
                        mutationMap.put(sampleId, features);
                    }


                    features.add(mut);
                }
            }


            return mutationMap;
        } catch (IOException e) {
            log.error("Error loading mutation file", e);
            throw new DataLoadException("IO Exception: " + e.toString(), locator.getPath());
        } finally {
            reader.close();
        }
    }

}
