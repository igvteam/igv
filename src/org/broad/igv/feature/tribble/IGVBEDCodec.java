/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.feature.tribble;

import org.broad.igv.feature.*;
import org.broad.tribble.util.ParsingUtils;

import java.util.LinkedHashMap;
import java.util.Map;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 20, 2009
 * Time: 10:15:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class IGVBEDCodec extends UCSCCodec {

    GFFParser.GFF3Helper tagHelper = new GFFParser.GFF3Helper();

    public BasicFeature decode(String nextLine) {

        if (nextLine.trim().length() == 0 || nextLine.startsWith("#") || nextLine.startsWith("track") ||
                nextLine.startsWith("browser")) {
            return null;
        }

        int tokenCount = ParsingUtils.splitWhitespace(nextLine, tokens);

        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.

        if (tokenCount < 2) {
            return null;
        }

        String chr = tokens[0];
        int start = Integer.parseInt(tokens[1]) - startBase;

        int end = start + 1;
        if (tokenCount > 2) {
            end = Integer.parseInt(tokens[2]) - startBase;
        }

        BasicFeature feature = new BasicFeature(chr, start, end);

        // The rest of the columns are optional.  Stop parsing upon encountering
        // a non-expected value

        // Name
        if (tokenCount > 3) {
            if (gffTags) {
                Map<String, String> atts = new LinkedHashMap();
                tagHelper.parseAttributes(tokens[3], atts);
                String name = tagHelper.getName(atts);
                //if (name == null) {
                //    name = tokens[3];
                //}
                feature.setName(name);

                String id = atts.get("ID");
                if (id != null) {
                    FeatureDB.put(id.toUpperCase(), feature);
                    feature.setIdentifier(id);
                } else {
                    feature.setIdentifier(name);
                }
                String alias = atts.get("Alias");
                if (alias != null) {
                    FeatureDB.put(alias.toUpperCase(), feature);
                }
                String geneSymbols = atts.get("Symbol");
                if (geneSymbols != null) {
                    String[] symbols = geneSymbols.split(",");
                    for (String sym : symbols) {
                        FeatureDB.put(sym.trim().toUpperCase(), feature);
                    }
                }

                String description = GFFParser.getDescription(atts);
                feature.setDescription(description);


            } else {
                String name = tokens[3].replaceAll("\"", "");
                feature.setName(name);
                feature.setIdentifier(name);
            }
        }

        // Score
        if (tokenCount > 4) {
            try {
                float score = Float.parseFloat(tokens[4]);
                feature.setScore(score);
            } catch (NumberFormatException numberFormatException) {

                // Unexpected, but does not invalidate the previous values.
                // Stop parsing the line here but keep the feature
                // Don't log, would just slow parsing down.
                return feature;
            }
        }

        // Strand
        if (tokenCount > 5) {
            String strandString = tokens[5].trim();
            char strand = (strandString.length() == 0)
                    ? ' ' : strandString.charAt(0);

            if (strand == '-') {
                feature.setStrand(Strand.NEGATIVE);
            } else if (strand == '+') {
                feature.setStrand(Strand.POSITIVE);
            } else {
                feature.setStrand(Strand.NONE);
            }
        }

        if (tokenCount > 8) {
            String colorString = tokens[8];
            feature.setColor(ParsingUtils.parseColor(colorString));
        }

        // Coding information is optional
        if (tokenCount > 11) {
            createExons(start, tokens, feature, chr, feature.getStrand());
        }

        return feature;
    }

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(String path) {
        return path.toLowerCase().endsWith(".bed");
    }


    private void createExons(int start, String[] tokens, BasicFeature gene, String chr,
                             Strand strand) throws NumberFormatException {

        int cdStart = Integer.parseInt(tokens[6]);
        int cdEnd = Integer.parseInt(tokens[7]);

        int exonCount = Integer.parseInt(tokens[9]);
        String[] exonSizes = new String[exonCount];
        String[] startsBuffer = new String[exonCount];
        ParsingUtils.split(tokens[10], exonSizes, ',');
        ParsingUtils.split(tokens[11], startsBuffer, ',');

        int exonNumber = (strand == Strand.NEGATIVE ? exonCount : 1);

        if (startsBuffer.length == exonSizes.length) {
            for (int i = 0; i < startsBuffer.length; i++) {
                int exonStart = start + Integer.parseInt(startsBuffer[i]);
                int exonEnd = exonStart + Integer.parseInt(exonSizes[i]);
                Exon exon = new Exon(chr, exonStart, exonEnd, strand);
                exon.setCodingStart(cdStart);
                exon.setCodingEnd(cdEnd);
                exon.setNumber(exonNumber);
                gene.addExon(exon);

                if (strand == Strand.NEGATIVE) {
                    exonNumber--;
                } else {
                    exonNumber++;
                }
            }
        }
    }

}