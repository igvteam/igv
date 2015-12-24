/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.feature.bionano;

import htsjdk.tribble.Feature;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * For the time being we bypass tribble and just use a parser
 */
public class SMAPParser {


    private SMAPParser() {

    }

    public static List<Feature> parseFeatures(ResourceLocator locator, Genome genome) throws IOException {

        BufferedReader br = null;

        ArrayList<Feature> features = new ArrayList<Feature>(1000);
        br = ParsingUtils.openBufferedReader(locator);
        String nextLine;
        while((nextLine = br.readLine()) != null) {

            if(nextLine.startsWith("#")) continue;

            String [] tokens = Globals.tabPattern.split(nextLine);

            String c1 = tokens[2];
            String c2 = tokens[3];
            int start = (int) Double.parseDouble(tokens[6]);
            int end = (int) Double.parseDouble(tokens[7]);

            if(c1.equals(c2)) {

                String chr = genome == null ? c1 : genome.getCanonicalChrName(c1);
                features.add(new SMAPFeature(chr, start, end, tokens));
            }
            else {
                String chr1 = genome == null ? c1 : genome.getCanonicalChrName(c1);
                features.add(new SMAPFeature(chr1, start, start + 1, tokens));

                String chr2 = genome == null ? c2 : genome.getCanonicalChrName(c2);
                features.add(new SMAPFeature(chr2, end, end + 1, tokens));
            }
        }

        return features;


    }


}
