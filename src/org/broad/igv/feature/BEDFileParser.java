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

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Parser for BED files
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/10
 */
public class BEDFileParser extends UCSCParser {

    private static Logger log = Logger.getLogger(BEDFileParser.class);

    private GFFParser.GFF3Helper tagHelper = new GFFParser.GFF3Helper();

    public BEDFileParser(Genome genome) {
        super(genome);
    }

    @Override
    /**
     * Test the contents of this file for the presence of features.  The first
     * 20 lines are parsed and if any features are found the method returns
     * true.  This is an imperfect test.  Other file type possibilities
     * should be eliminated before using this method
     *
     *
     * @param locator
     * @return true if a feature file
     */
    public boolean isFeatureFile(ResourceLocator locator) {
        String ext = getStrippedFilename(locator.getPath());
        return ext.endsWith("bed");

    }

    public BasicFeature parseLine(String[] tokens, int tokenCount) {


        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.

        if (tokenCount < 2) {
            return null;
        }

        String chr = genome == null ? tokens[0] : genome.getChromosomeAlias(tokens[0]);
        int start = Integer.parseInt(tokens[1]) - startBase;

        int end = start + 1;
        if (tokenCount > 2) {
            end = Integer.parseInt(tokens[2]);
        }

        BasicFeature feature = null;

        boolean isSpliceJunction = false;
        //dhmay adding.  This is a hacky way to determine whether the BED file we're reading is a
        //splice junction BED file.
        //todo: some refactoring that allows this hack to be removed
        if (trackProperties != null) {
            Class rendererClass = trackProperties.getRendererClass();
            if (rendererClass != null && rendererClass.isAssignableFrom(SpliceJunctionRenderer.class))
                isSpliceJunction = true;
        }

        if (isSpliceJunction)
            feature = new SpliceJunctionFeature(chr, start, end);
        else
            feature = new BasicFeature(chr, start, end);

        // The rest of the columns are optional.  Stop parsing upon encountering
        // a non-expected value

        // Name
        if (tokenCount > 3) {
            if (gffTags) {
                Map<String, String> atts = new HashMap();
                tagHelper.parseAttributes(tokens[3], atts);
                String name = tagHelper.getName(atts);
                if (name == null) {
                    name = tokens[3];
                }
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
                if(geneSymbols != null) {
                    String [] symbols = geneSymbols.split(",");
                    for(String sym : symbols) {
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
                //todo: some refactoring that allows this hack to be removed
                if (isSpliceJunction)
                    ((SpliceJunctionFeature) feature).setJunctionDepth((int) score);
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

        // thick start & end
        if (tokenCount > 7) {
            feature.setThickStart(Integer.parseInt(tokens[6]));
            feature.setThickEnd(Integer.parseInt(tokens[7]));
        }

        if (tokenCount > 8) {
            String[] rgb = new String[3];
            int nTokens = ParsingUtils.split(tokens[8].replaceAll("\"", ""), rgb, ',');
            feature.setColor(rgb, nTokens);
        }

        // Coding information is optional, except in the case of splice junctions
        if (tokenCount > 11) {
            createExons(start, tokens, feature, chr, feature.getStrand());

            //todo: some refactoring that allows this hack to be removed
            if (isSpliceJunction) {
                SpliceJunctionFeature junctionFeature = (SpliceJunctionFeature) feature;

                List<Exon> exons = feature.getExons();

                junctionFeature.setJunctionStart(start + exons.get(0).getLength());
                junctionFeature.setJunctionEnd(end - exons.get(1).getLength());

                new StringBuilder();
            }
        }


        return feature;
    }

    private void createExons(int start, String[] tokens, BasicFeature gene, String chr,
                             Strand strand) throws NumberFormatException {

        int cdStart = Integer.parseInt(tokens[6]) - startBase;
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

        // Do not compute reading shifts for bed files.  Just assume they are non-coding.
    }
}
