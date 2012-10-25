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
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.MutationTrack;
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

    public static boolean isMutationAnnotationFile(ResourceLocator locator) throws IOException {
        return MUTCodec.isMutationAnnotationFile(locator.getPath());

    }

    public List<FeatureTrack> loadMutationTracks(ResourceLocator locator, Genome genome) {
        this.locator = locator;
        this.genome = genome;

        List<FeatureTrack> tracks = new ArrayList();
        Map<String, List<org.broad.tribble.Feature>> features = loadMutations();
        for (String sampleId : features.keySet()) {
            String id = locator.getPath() + "_" + sampleId;
            MutationTrack track = new MutationTrack(locator, id, new FeatureCollectionSource(features.get(sampleId), genome));
            track.setSquishedRowHeight(5);
            track.setExpandedRowHeight(15);
            track.setHeight(15);
            track.setName(sampleId);

            // Override default minimum height (10 for feature tracks).
            track.setMinimumHeight(0);
            tracks.add(track);
        }
        //Just to make sure we have no memory
        this.locator = null;
        this.genome = null;
        return tracks;
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

            MUTCodec codec = new MUTCodec(locator.getPath(), genome);

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
