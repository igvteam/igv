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

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class RNAIGeneScoreParser {

    public enum Type {
        GENE_SCORE, POOLED
    }

    private static Logger log = Logger.getLogger(RNAIGeneScoreParser.class);
    private String filename;

    private int maxColumn = -1;
    private int batchColumn = -1;
    private int conditionColumn = -1;
    private int geneColumn = -1;
    private int hairpinColumn = -1;
    private int scoreColumn = -1;
    private int confidenceColumn = -1;
    private Genome genome;

    public RNAIGeneScoreParser(String filename, Type type, Genome genome) {

        this.genome = genome;
        this.filename = filename;
        if (type == Type.GENE_SCORE) {
            batchColumn = 0;
            conditionColumn = 2;
            geneColumn = 5;
            hairpinColumn = 12;
            scoreColumn = 8;
            confidenceColumn = 9;
            maxColumn = 12;
        } else {
            batchColumn = -1;
            conditionColumn = -1;
            geneColumn = 0;
            hairpinColumn = 1;
            scoreColumn = 3;
            confidenceColumn = 2;
            maxColumn = 3;
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public Collection<RNAIDataSource> parse() {

        BufferedReader reader = null;

        try {
            log.debug("Loading data for: " + filename);
            reader = new BufferedReader(new FileReader(filename));

            // Parse comments
            parseAttributes(reader);

            // Parse header
            parseHeaderRow(reader);


            // Parse data
            String nextLine = reader.readLine();

            // Skip empty and comment lines
            while ((nextLine = reader.readLine()) != null
                    && (nextLine.startsWith("#") || (nextLine.length() == 0))) {
            }

            Map<String, RNAIDataSource> dataSources = new HashMap();
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                if (tokens.length > maxColumn) {
                    try {
                        String batchId = (batchColumn < 0) ? "" : tokens[batchColumn].trim();
                        String geneName = tokens[geneColumn].trim().toUpperCase();
                        NamedFeature gene = FeatureDB.getFeature(geneName);

                        if (gene != null) {
                            float geneScore = Float.NaN;
                            try {
                                geneScore = Float.parseFloat(tokens[scoreColumn]);

                            }
                            catch (NumberFormatException numberFormatException) {

                                // Nothing to do -- expected condition.  Indicates no score for this
                                // gene
                            }

                            int confidence = 0;
                            try {
                                confidence = Integer.parseInt(tokens[confidenceColumn]);
                            }
                            catch (NumberFormatException numberFormatException) {

                                // Nothing to do -- expected condition.  This will occur when the
                                // score is blank.
                            }

                            int numberOfHairpins = 0;
                            try {
                                numberOfHairpins = Integer.parseInt(tokens[hairpinColumn]);
                            }
                            catch (NumberFormatException numberFormatException) {

                                // Nothing to do -- expected condition.  This will coincide with
                                // a blank gene score
                            }

                            // Make batch_conditon key
                            // Rules from Jessee -- ignore conditions starting with *.  None, standard,
                            // and blank are all equivalent.
                            String cond = (conditionColumn < 0)
                                    ? "" : tokens[conditionColumn].trim();
                            if (!cond.startsWith("*")) {
                                if (cond.equals("None") || cond.equals("Standard")) {
                                    cond = "";
                                }

                                String batchCond = batchId + "_" + cond;

                                RNAIDataSource ds = dataSources.get(batchCond);

                                // List<RNAIGeneScore> dataPoints = dataPointsMap.get(batchCond);

                                if (ds == null) {
                                    ds = new RNAIDataSource(batchId, cond, genome);
                                    dataSources.put(batchCond, ds);
                                }
                                ds.addGeneScore(new RNAIGeneScore(batchCond, gene, geneScore,
                                        confidence, numberOfHairpins));
                            }

                        } else {

                            // todo -- handle unknown gene
                            log.info("Unknown gene: " + geneName);
                        }
                    }
                    catch (Exception ex) {
                        log.error("Skipping line: " + nextLine, ex);
                    }
                }
            }

            return dataSources.values();

        }
        catch (FileNotFoundException e) {
            log.error("RNAI file not found: " + filename);
            throw new RuntimeException("File not found");
        }
        catch (IOException e) {
            log.error(filename, e);
            throw new RuntimeException("Error parsing file.", e);
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();

                }
                catch (IOException iOException) {
                }
            }
        }
    }

    /**
     * Parse the attributes from the comment section and annotate the
     */
    private void parseAttributes(BufferedReader reader) throws IOException {

        // TODO -- parse comments
    }

    private void parseHeaderRow(BufferedReader reader) throws IOException {

        // Nothing to do here.  Column positions are fixed.
    }
}
