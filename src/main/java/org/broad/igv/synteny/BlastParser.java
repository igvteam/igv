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

package org.broad.igv.synteny;

//~--- JDK imports ------------------------------------------------------------

import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.Strand;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class BlastParser {

    /**
     * Method description
     * <p/>
     * queryID,subjectID,percentIden,alignmentLength,
     * mismatches,gapOpenings,queryStart,queryStop,subjectStart,subjectStop,eval,bitScore
     * <p/>
     * 1	9	87.73	864	106	0	100842	101705	222166	221303	0	872
     *
     * @param file
     * @return
     */
    public List<BlastMapping> parse(String file) {

        List<BlastMapping> mappings = new ArrayList();

        BufferedReader reader = null;

        boolean tabblastn = file.contains(".tabblastn");

        int queryChrIdx = 0;
        int queryStartIdx = tabblastn ? 6 : 1;
        int queryEndIdx = tabblastn ? 7 : 2;
        int queryStrandIdx = tabblastn ? -1 : 3;

        int subjectChrIdx = tabblastn ? 1 : 4;
        int subjectStartIdx = tabblastn ? 8 : 5;
        int subjectEndIdx = tabblastn ? 9 : 6;
        int subjectStrandIdx = tabblastn ? -1 : 7;

        int nTokens = tabblastn ? 10 : 8;

        long lineCount = 0;
        String nextLine = null;
        boolean dataParsed = false;

        try {
            reader = new BufferedReader(new FileReader(file));

            while ((nextLine = reader.readLine()) != null) {
                lineCount++;
                String[] tokens = nextLine.split("\t");
                if (tokens.length >= nTokens) {
                    dataParsed = true;
                    String queryCht = tokens[queryChrIdx];

                    int queryStart;
                    int queryEnd;
                    try {
                        queryStart = Integer.parseInt(tokens[queryStartIdx]);
                        queryEnd = Integer.parseInt(tokens[queryEndIdx]);
                    }
                    catch (NumberFormatException ne) {
                        throw new RuntimeException("Non-numeric value found in either start or end column");
                    }
                    Strand queryStrand = Strand.NONE;
                    if (!tabblastn) {
                        int str = Integer.parseInt(tokens[queryStrandIdx]);
                        queryStrand = str == 1 ? Strand.POSITIVE : Strand.NEGATIVE;
                    }

                    String subjectChr = tokens[subjectChrIdx];

                    int subjectStart;
                    int subjectEnd;
                    try {
                        subjectStart = Integer.parseInt(tokens[subjectStartIdx]);
                        subjectEnd = Integer.parseInt(tokens[subjectEndIdx]);
                    }
                    catch (NumberFormatException ne) {
                        throw new RuntimeException("Non-numeric value found in " +
                                " either subject start or end column");
                    }
                    Strand subjectStrand = Strand.NONE;
                    if (!tabblastn) {
                        try {
                            int str = Integer.parseInt(tokens[subjectStrandIdx]);
                            subjectStrand = str == 1 ? Strand.POSITIVE : Strand.NEGATIVE;
                        }
                        catch (NumberFormatException ne) {
                            throw new RuntimeException("Non-numeric value found in " +
                                    " either subject strand column");
                        }
                    }

                    float percentIden = tabblastn ? Float.parseFloat(tokens[2]) : 100.0f;

                    BlastMapping.Block queryBlock = new BlastMapping.Block(queryCht, queryStart, queryEnd, queryStrand);
                    BlastMapping.Block subjectBlock = new BlastMapping.Block(subjectChr, subjectStart, subjectEnd, subjectStrand);

                    mappings.add(new BlastMapping(queryBlock, subjectBlock, percentIden));
                }
            }

            if (!dataParsed) {
                throw new RuntimeException("Data not loaded. Number of columns is less than " + nTokens + ".");
            }

        } catch (Exception ex) {
            if (lineCount != 0) {
                throw new ParserException(ex.getMessage(), lineCount, nextLine);
            } else {
                throw new RuntimeException(ex);
            }
        } finally {
            if (reader != null) {
                try {
                    reader.close();

                } catch (IOException iOException) {
                }
            }
        }

        return mappings;
    }
}
