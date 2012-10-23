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
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.MutationTrack;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.MultiMap;
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
public class MutationParser {

    private static Logger log = Logger.getLogger(MutationParser.class);
    private int chrColumn;
    private int startColumn;
    private int endColumn;
    private int sampleColumn;
    private int typeColumn;
    private int refAlleleColumn;
    private int tumorAllele1Column;
    private int tumorAllele2Column;

    private boolean isMAF;
    private String[] headers = null;

    private ResourceLocator locator = null;
    private Genome genome = null;

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

    private void readHeader(String[] headers) {
        isMAF = headers.length > 15 && headers[0].equalsIgnoreCase("Hugo_Symbol");
        setColumns(isMAF);
        this.headers = headers;
    }

    private Mutation decode(String[] tokens) {
        String chr = genome.getChromosomeAlias(tokens[chrColumn].trim());

        int start;
        try {
            start = Integer.parseInt(tokens[startColumn].trim());
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (startColumn + 1) + " must be a numeric value.", locator.getPath());
        }

        int end;
        try {
            end = Integer.parseInt(tokens[endColumn].trim());
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (endColumn + 1) + " must be a numeric value.", locator.getPath());
        }


        // MAF files use the 1-based inclusive convention for coordinates.  The convention is not
        // specified for MUT files, and it appears both conventions have been used.  We can detect
        // the convention used for single base mutations by testing start == end.
        if (isMAF || (start == end)) {
            start--;
        }

        String sampleId = tokens[sampleColumn].trim();
        String type = tokens[typeColumn].trim();

        MultiMap<String, String> attributes = new MultiMap();
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

        if (refAlleleColumn > 0) {
            mut.setRefAllele(tokens[refAlleleColumn].trim());
        }
        if (tumorAllele1Column > 0) {
            mut.setAltAllele1(tokens[tumorAllele1Column].trim());
        }
        if (tumorAllele2Column > 0) {
            mut.setAltAllele2(tokens[tumorAllele2Column].trim());
        }

        return mut;
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

            Map<String, List<org.broad.tribble.Feature>> mutationMap = new LinkedHashMap();

            reader = ParsingUtils.openAsciiReader(locator);

            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("#")) continue;

                String[] tokens = nextLine.split("\t");
                if (tokens.length > 4) {

                    if (this.headers == null) {
                        readHeader(tokens);
                        continue;
                    }

                    Mutation mut = decode(tokens);

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

    private void setColumns(boolean isMAF) {

        if (isMAF) {
            chrColumn = 4;
            startColumn = 5;
            endColumn = 6;
            sampleColumn = 15;
            typeColumn = 8;
            refAlleleColumn = 10;
            tumorAllele1Column = 11;
            tumorAllele2Column = 12;
        } else {
            chrColumn = 0;
            startColumn = 1;
            endColumn = 2;
            sampleColumn = 3;
            typeColumn = 4;
            refAlleleColumn = -1;
            tumorAllele1Column = -1;
            tumorAllele2Column = -1;
        }
    }
}
