/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.broad.igv.feature.genome.Genome;

/**
 * Parser for the UCSC snp tracks (e.g. snp128.txt).
 *
 * @author jrobinso
 */
public class SnpFileParser extends AbstractFeatureParser {

    public SnpFileParser(Genome genome) {
        super(genome);
    }

    protected IGVFeature parseLine(String line) throws NumberFormatException {

        String[] tokens = line.replaceAll("\"", "").split("\t");

        int tokenCount = tokens.length;
        int i = 0;
        String chr = tokens[i++];
        int start = Integer.parseInt(tokens[i++]);
        int end = Integer.parseInt(tokens[i++]);
        BasicFeature feature = new BasicFeature(chr, start, end);
        if (tokenCount > 3) {
            String name = tokens[i++].replaceAll("\"", "");
            feature.setName(name);
        }
        if (tokenCount > 4) {
            float score = Float.parseFloat(tokens[i++]);
            feature.setScore(score);
        }
        if (tokenCount > 5) {
            String strandString = tokens[i++].trim();
            char strand = strandString.length() == 0 ? ' ' : strandString.charAt(0);

            if (strand == '-') {
                feature.setStrand(Strand.NEGATIVE);
            } else if (strand == '+') {
                feature.setStrand(Strand.POSITIVE);
            } else {
                feature.setStrand(Strand.NONE);
            }
        }

        // TODO -- rest of file

        return feature;
    }

}
