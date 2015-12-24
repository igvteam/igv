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
////SmapEntryID	QryContigID	RefcontigID1	RefcontigID2	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Confidence	Type	XmapID1	XmapID2	LinkID	QryStartIdx	QryEndIdx	RefStartIdx	RefEndIdx

    public static List<Feature> parseFeatures(ResourceLocator locator, Genome genome) throws IOException {

        BufferedReader br = null;

        ArrayList<Feature> features = new ArrayList<Feature>(1000);
        br = ParsingUtils.openBufferedReader(locator);
        String nextLine;
        String[] headers = null;
        int refContig1 = -1, refContig2 = -1, refStart = -1, refEnd = -1, confidence = -1, type = -1;



        while ((nextLine = br.readLine()) != null) {

            if (nextLine.startsWith("#h")) {
                headers = Globals.tabPattern.split(nextLine.substring(3));
                for (int i = 0; i < headers.length; i++) {
                    String h = headers[i];
                    if (h.equals("RefcontigID1")) {
                        refContig1 = i;
                    } else if (h.equals("RefcontigID2")) {
                        refContig2 = i;
                    } else if (h.equals("RefStartPos")) {
                        refStart = i;
                    } else if (h.equals("RefEndPos")) {
                        refEnd = i;
                    } else if (h.equals("Confidence")) {
                        confidence = i;
                    } else if (h.equals("Type")) {
                        type = i;
                    }
                }
            } else if (nextLine.startsWith("#")) {
                continue;
            } else {
                if (headers == null) {
                    throw new RuntimeException("Never saw #h line");
                }

                String[] tokens = Globals.tabPattern.split(nextLine);

                String c1 = tokens[refContig1];
                String c2 = tokens[refContig2];
                int start = (int) Double.parseDouble(tokens[refStart]);
                int end = (int) Double.parseDouble(tokens[refEnd]);
                double conf = Double.parseDouble(tokens[confidence]);
                String t = tokens[type];


                if (c1.equals(c2)) {

                    String chr = genome == null ? c1 : genome.getCanonicalChrName(c1);
                    features.add(new SMAPFeature(chr, start, end, conf, t, headers, tokens));
                } else {
                    String chr1 = genome == null ? c1 : genome.getCanonicalChrName(c1);
                    features.add(new SMAPFeature(chr1, start, end, conf, t, headers, tokens));

                    String chr2 = genome == null ? c2 : genome.getCanonicalChrName(c2);
                    features.add(new SMAPFeature(chr2, start, end, conf, t, headers, tokens));
                }
            }
        }

        return features;


    }


}
