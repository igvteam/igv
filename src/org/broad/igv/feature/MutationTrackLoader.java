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
package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.MUTCodec;
import org.broad.igv.track.*;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Parses a ".mut" file
 *
 * @author jrobinso
 */
public class MutationTrackLoader {

    private static Logger log = Logger.getLogger(MutationTrackLoader.class);
    private ResourceLocator locator = null;
    private Genome genome = null;
    MUTCodec codec;

    public static boolean isMutationAnnotationFile(ResourceLocator locator) throws IOException {
        return MUTCodec.isMutationAnnotationFile(locator.getPath());
    }

    public List<FeatureTrack> loadMutationTracks(ResourceLocator locator, Genome genome) throws IOException {

        this.locator = locator;
        this.genome = genome;

        boolean indexed = isIndexed(locator.getPath(), genome);

        List<FeatureTrack> tracks = new ArrayList();

        if (indexed) {
            String[] samples = getCodec().getSamples();
            MutationDataManager dataManager = new MutationDataManager(locator.getPath(), genome);
            for (String sampleId : samples) {
                String id = locator.getPath() + "_" + sampleId;
                FeatureSource<Mutation> featureSource = new MutationFeatureSource(sampleId, dataManager);
                MutationTrack track = new MutationTrack(locator, id, featureSource);
                tracks.add(track);
                track.setName(sampleId);
            }

        } else {
            Map<String, List<org.broad.tribble.Feature>> features = loadMutations();
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
     * Test to see if a usable index exists.  In addition to the index, mutaion files have an additional requirement
     * that samples be specified in a header directive.
     *
     * @param path
     * @return
     */
    private boolean isIndexed(String path, Genome genome) {
        if (!TrackLoader.isIndexed(path, genome)) return false;

        try {
            String[] samples = getCodec().getSamples();
            return samples != null && samples.length > 0;
        } catch (Exception e) {
            log.error("Error creating codec for: " + path, e);
            return false;
        }

    }

    /**
     * Return a map of runId -> list of mutation objects.   The "runId" field
     * is the track identifier (name) for mutation files.
     *
     * @return
     */
    private Map<String, List<org.broad.tribble.Feature>> loadMutations() {
        AsciiLineReader reader = null;
        String nextLine = null;

        try {

            if (codec == null) codec = new MUTCodec(locator.getPath(), genome);

            Map<String, List<org.broad.tribble.Feature>> mutationMap = new LinkedHashMap();

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
                    List<org.broad.tribble.Feature> features = mutationMap.get(sampleId);
                    if (features == null) {
                        features = new ArrayList();
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
