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

package org.broad.igv.data.rnai;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.readers.AsciiLineReader;

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
