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

package org.broad.igv.data.rnai;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.data.expression.ProbeToLocusMap;
import org.broad.igv.exceptions.LoadResourceFromServerException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: nazaire
 * Date: Feb 25, 2009
 */
public class RNAIGCTDatasetParser {

    private static Logger log = Logger.getLogger(RNAIGCTDatasetParser.class);
    private ResourceLocator dataFileLocator;
    private int dataStartColumn = 2;
    private int descriptionColumn = 1;
    Genome genome;


    /**
     * Constructs ...
     *
     * @param gctFile
     */
    public RNAIGCTDatasetParser(ResourceLocator gctFile, Genome genome) {
        this.dataFileLocator = gctFile;
        this.genome = genome;
        dataStartColumn = 2;
    }

    public Collection<RNAIDataSource> parse() {
        // Create a buffer for the string split utility.  We use  a custom utility as opposed
        AsciiLineReader reader = null;
        List dataSources = null;
        String nextLine = null;
        InputStream probeMappingStream = null;

        try {
            reader = ParsingUtils.openAsciiReader(dataFileLocator);

            String headerLine = null;

            // Skip header rows

            nextLine = reader.readLine();
            nextLine = reader.readLine();

            headerLine = reader.readLine();

            // Parse column headings
            int skip = 1;

            String[] tokens = Globals.tabPattern.split(headerLine, -1);
            int nTokens = tokens.length;

            String description = (nTokens > descriptionColumn)
                    ? new String(tokens[descriptionColumn]) : null;

            int nColumns = (nTokens - dataStartColumn) / skip;
            String[] columnHeadings = new String[nColumns];
            for (int i = 0; i < nColumns; i++) {
                String heading = tokens[dataStartColumn + i * skip].replace('\"', ' ').trim();
                columnHeadings[i] = heading;
            }


            Map<String, String[]> rnaiProbeMap = getProbeMap();

            HashMap<String, HashMap<String, Float>> sampleGeneScoreMap = new HashMap();
            while ((nextLine = reader.readLine()) != null) {
                tokens = Globals.tabPattern.split(nextLine, -1);
                nTokens = tokens.length;
                String probeId = new String(tokens[0]);
                float[] values = new float[nColumns];

                String[] identifiers = (String[]) rnaiProbeMap.get(probeId);
                String identifier = null;
                if (identifiers == null || identifiers.length == 0) {
                    log.info("Could not find mapping for: " + probeId);
                    continue;
                } else {
                    identifier = identifiers[0];
                }

                NamedFeature gene = FeatureDB.getFeature(identifier.toUpperCase());
                if (gene == null) {
                    log.debug("Unknown identifier: " + identifier);
                    continue;
                }

                for (int i = 0; i < nColumns; i++) {
                    try {
                        int dataIndex = dataStartColumn + i * skip;

                        // If we are out of value tokens, or the cell is blank, assign NAN to the cell.
                        if ((dataIndex >= nTokens) || (tokens[dataIndex].length() == 0)) {
                            values[i] = Float.NaN;
                        } else {
                            values[i] = Float.parseFloat(tokens[dataIndex]);
                        }

                        String sample = columnHeadings[i];
                        RNAIHairpinValue hairpin = new RNAIHairpinValue(probeId, values[i]);
                        RNAIHairpinCache.getInstance().addHairpinScore(sample, gene.getName(),
                                hairpin);

                        HashMap<String, Float> geneScoreMap = sampleGeneScoreMap.get(sample);

                        if (geneScoreMap == null) {
                            geneScoreMap = new HashMap();
                            sampleGeneScoreMap.put(sample, geneScoreMap);
                        }

                        Float geneScore = geneScoreMap.get(gene.getName());
                        if (geneScore == null) {
                            geneScore = values[i];
                            geneScoreMap.put(gene.getName(), geneScore);
                        } else {

                            geneScore = new Float(Math.min(values[i], geneScore.floatValue()));
                            geneScoreMap.put(gene.getName(), geneScore);
                        }
                    } catch (NumberFormatException numberFormatException) {

                        // This is an expected condition.  IGV uses NaN to
                        // indicate non numbers (missing data values)
                        values[i] = Float.NaN;
                    }
                }
            }

            dataSources = computeGeneScores(sampleGeneScoreMap);

        } catch (IOException ex) {
            log.error("Error parsing RNAi file", ex);
            throw new RuntimeException(ex);
        } finally {
            if (probeMappingStream != null) {
                try {
                    probeMappingStream.close();
                } catch (IOException e) {
                    log.error("Error closing probe mapping stream", e);
                }
            }
            if (reader != null) {
                reader.close();
            }
        }

        return dataSources;
    }

    private List computeGeneScores(HashMap<String, HashMap<String, Float>> sampleGeneScoreMap) {
        int confidence = 3;
        List dataSources = new ArrayList();
        Iterator samplesIt = sampleGeneScoreMap.keySet().iterator();
        while (samplesIt.hasNext()) {
            String sample = (String) samplesIt.next();
            HashMap geneMap = sampleGeneScoreMap.get(sample);
            RNAIDataSource ds = new RNAIDataSource(sample, "", genome);
            Iterator geneScoreIt = geneMap.keySet().iterator();
            while (geneScoreIt.hasNext()) {
                String gene = (String) geneScoreIt.next();
                Float score = (Float) geneMap.get(gene);
                int numHairpins;
                Collection hairpins = RNAIHairpinCache.getInstance().getHairpinScores(sample, gene);
                if (hairpins == null) {
                    numHairpins = 0;
                } else {
                    numHairpins = hairpins.size();
                }
                ds.addGeneScore(new RNAIGeneScore(sample,
                        FeatureDB.getFeature(gene), score.floatValue(), numHairpins));
            }
            dataSources.add(ds);
        }
        return dataSources;
    }


    private final static String RNAI_MAPPING_FILE = "http://www.broadinstitute.org/igv/resources/probes/rnai/RNAI_probe_mapping.txt.gz";
    private static Map<String, String[]> rnaiProbeMap = null;

    private synchronized static Map<String, String[]> getProbeMap() throws IOException {
        if (rnaiProbeMap == null) {
            rnaiProbeMap = Collections.synchronizedMap(new HashMap<String, String[]>(20000));
            URL url = new URL(RNAI_MAPPING_FILE);

            InputStream probeMappingStream = null;
            try {
                probeMappingStream = new GZIPInputStream(HttpUtils.getInstance().openConnectionStream(url));
                AsciiLineReader br = new AsciiLineReader(probeMappingStream);

                ProbeToLocusMap.getInstance().loadMapping(br, rnaiProbeMap);
            } catch (Exception e) {
                throw new LoadResourceFromServerException(e.getMessage(), RNAI_MAPPING_FILE, e.getClass().getSimpleName());
            } finally {
                if (probeMappingStream != null) {
                    probeMappingStream.close();
                }
            }
        }

        return rnaiProbeMap;
    }
}
