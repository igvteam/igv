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

package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 4/16/13
 *         Time: 2:25 PM
 */
public class ChromSizesParser {

    private static Logger log = Logger.getLogger(ChromSizesParser.class);

    public static List<Chromosome> parse(String path) throws IOException {

        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(path);
            List<Chromosome> chromosomes = new ArrayList<Chromosome>();
            String nextLine;
            int idx = 0;
            while ((nextLine = br.readLine()) != null) {

                String[] tokens = Globals.whitespacePattern.split(nextLine);
                if (tokens.length >= 2) {
                    String chr = tokens[0];
                    int length = Integer.parseInt(tokens[1]);
                    chromosomes.add(new Chromosome(idx, chr, length));
                    idx++;
                } else {
                    log.info("Unexpected # of tokens at line: " + nextLine);
                }

            }
            return chromosomes;
        } finally {
            if (br != null) br.close();
        }

    }


}
