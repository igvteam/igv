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
import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * @author jrobinso
 */
public class RNAIHairpinParser {

    private static Logger log = Logger.getLogger(RNAIHairpinParser.class);
    private String filename;

    /**
     * Constructs ...
     *
     * @param filename
     */
    public RNAIHairpinParser(String filename) {

        this.filename = filename;
    }

    /**
     * Method description
     */
    public void parse() {


        AsciiLineReader reader = null;

        try {
            log.debug("Loading data for: " + filename);
            reader = new AsciiLineReader(new FileInputStream(filename));

            // Parse comments
            parseAttributes(reader);

            // Parse header
            parseHeaderRow(reader);

            // Parse data
            String nextLine = reader.readLine();

            while ((nextLine = reader.readLine()) != null) {
                if (!nextLine.startsWith("#")) {

                    String[] tokens = Globals.tabPattern.split(nextLine, -1);
                    int nTokens = tokens.length;
                    if (nTokens > 11) {
                        try {
                            String batchId = new String(tokens[0].trim());

                            String hairpinName = new String(tokens[4].trim().toUpperCase());

                            float scoreMean = Float.NaN;
                            try {
                                scoreMean = Float.parseFloat(tokens[8]);

                            } catch (NumberFormatException numberFormatException) {

                                // Nothing to do -- expected condition
                            }

                            float scoreSTD = 0;
                            try {
                                scoreSTD = Float.parseFloat(tokens[9]);
                            } catch (NumberFormatException numberFormatException) {

                                // Nothing to do -- expected condition
                            }

                            String gene = new String(tokens[13].trim().toUpperCase());

                            RNAIHairpinValue score = new RNAIHairpinValue(hairpinName, scoreMean,
                                    scoreSTD);

                            // Make batch_conditon key
                            // Rules from Jessee -- ignore conditions starting with *.  None, standard,
                            // and blank are all equivalent.
                            String cond = tokens[2].trim();
                            if (!cond.startsWith("*")) {
                                if (cond.equals("None") || cond.equals("Standard")) {
                                    cond = "";
                                }

                                String batchCond = batchId + "_" + cond;

                                RNAIHairpinCache.getInstance().addHairpinScore(batchCond, gene,
                                        score);
                            }
                        } catch (Exception ex) {
                            log.error("Skipping line: " + nextLine, ex);
                        }
                    }
                }
            }

        } catch (FileNotFoundException e) {
            log.error("RNAI file not found: " + filename);
            throw new RuntimeException("File not found: " + filename);

        } catch (IOException e) {
            log.error(filename, e);
            throw new RuntimeException("Error parsing file " + filename, e);
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    /**
     * Parse the attributes from the comment section and annotate the
     */
    private void parseAttributes(AsciiLineReader reader) throws IOException {

        // TODO -- parse comments
    }

    private void parseHeaderRow(AsciiLineReader reader) throws IOException {

        // Nothing to do here.  Column positions are fixed.
    }
}
