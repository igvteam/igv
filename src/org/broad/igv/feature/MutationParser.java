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
package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.MutationTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

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
public class MutationParser {

    private static Logger log = Logger.getLogger(MutationParser.class);
    private int chrColumn;
    private int startColumn;
    private int endColumn;
    private int sampleColumn;
    private int typeColumn;
    private int refAlleleColumn;
    private int tumorAlleleColumn;

    public static boolean isMutationAnnotationFile(ResourceLocator locator) throws IOException {
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(locator);
            if (reader == null) {
                return false;
            }

            String nextLine;
            while ((nextLine = reader.readLine()) != null && nextLine.startsWith("#")) {
                if (nextLine.startsWith("#version")) {
                    return true;
                }
            }
            if (nextLine == null) return false;

            String[] tokens = nextLine.split("\t");
            return tokens.length > 15 && tokens[0].equalsIgnoreCase("Hugo_Symbol");
        } catch (Exception e) {
            log.error("", e);
            throw new DataLoadException("Error checking for MAF file type: " + e.toString(), locator.getPath());
        } finally {
            if (reader != null) {
                reader.close();
            }
        }


    }

    public List<FeatureTrack> loadMutationTracks(ResourceLocator locator, Genome genome) {

        List<FeatureTrack> tracks = new ArrayList();
        Map<String, List<org.broad.tribble.Feature>> features = loadMutations(locator, genome);
        for (String sampleId : features.keySet()) {
            String id = locator.getPath() + "_" + sampleId;
            MutationTrack track = new MutationTrack(locator, id, new FeatureCollectionSource(features.get(sampleId), genome));
            track.setSquishedRowHeight(5);
            track.setExpandedRowHeight(15);
            track.setHeight(15);
            track.setName(sampleId);

            // Overrid default minimum height (10 for feature tracks).
            track.setMinimumHeight(0);
            tracks.add(track);
        }
        return tracks;
    }

    /**
     * Return a map of runId -> list of mutation objects.   The "runId" field
     * is the track identifier (name) for mutation files.
     *
     * @param locator
     * @return
     */
    private Map<String, List<org.broad.tribble.Feature>> loadMutations(ResourceLocator locator, Genome genome) {
        AsciiLineReader reader = null;
        String nextLine = null;

        try {

            Map<String, List<org.broad.tribble.Feature>> mutationMap = new LinkedHashMap();

            reader = ParsingUtils.openAsciiReader(locator);

            String[] headers = null;
            int lineNumber = 1;
            boolean isMAF = false;
            while ((nextLine = reader.readLine()) != null) {
                lineNumber++;

                if(nextLine.startsWith("#")) continue;

                String[] tokens = nextLine.split("\t");
                if (tokens.length > 4) {

                    if (headers == null) {
                        headers = tokens;
                        isMAF = headers.length > 15 && headers[0].equalsIgnoreCase("Hugo_Symbol");
                        setColumns(headers, isMAF);
                        continue;
                    }


                    String chr = genome.getChromosomeAlias(tokens[chrColumn].trim());

                    int start;
                    try {
                        start = Integer.parseInt(tokens[startColumn].trim());
                        if (isMAF) {
                            start--;
                        }
                    } catch (NumberFormatException e) {
                        throw new DataLoadException("Column " + (startColumn + 1) + " must be a numeric value.", locator.getPath());
                    }

                    int end;
                    try {
                        end = Integer.parseInt(tokens[endColumn].trim());
                    } catch (NumberFormatException e) {
                        throw new DataLoadException("Column " + (endColumn + 1) + " must be a numeric value.", locator.getPath());
                    }

                    String sampleId = tokens[sampleColumn].trim();
                    String type = tokens[typeColumn];

                    LinkedHashMap<String, String> attributes = new LinkedHashMap();
                    int n = Math.min(headers.length, tokens.length);
                    for (int i = 0; i < n; i++) {
                        String key = headers[i];
                        String value = tokens[i];
                        if (value.length() > 0) {
                            attributes.put(key, value);
                        }
                    }



                    Mutation mut = new Mutation(sampleId, chr, start, end, type);
                    mut.setAttributes(attributes);

                    if(refAlleleColumn > 0) {
                        mut.setRefAllele(tokens[refAlleleColumn].trim());
                    }
                    if(tumorAlleleColumn > 0) {
                        mut.setAltAllele(tokens[tumorAlleleColumn].trim());
                    }

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
            log.error("", e);
            throw new DataLoadException("IO Exception: " + e.toString(), locator.getPath());
        } finally {
            reader.close();
        }
    }

    private void setColumns(String[] headings, boolean isMAF) {

        if (isMAF) {
            chrColumn = 4;
            startColumn = 5;
            endColumn = 6;
            sampleColumn = 15;
            typeColumn = 8;
            refAlleleColumn = 10;
            tumorAlleleColumn = 11;
        } else {
            chrColumn = 0;
            startColumn = 1;
            endColumn = 2;
            sampleColumn = 3;
            typeColumn = 4;
            refAlleleColumn = -1;
            tumorAlleleColumn = -1;
        }
    }
}
