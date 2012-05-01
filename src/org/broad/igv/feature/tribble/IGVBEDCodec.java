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

import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.tribble.util.ParsingUtils;

import java.net.URLEncoder;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 20, 2009
 * Time: 10:15:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class IGVBEDCodec extends UCSCCodec {

    static final Pattern BR_PATTERN = Pattern.compile("<br>");
    static final Pattern EQ_PATTERN = Pattern.compile("=");

    Genome genome;

    public IGVBEDCodec() {
        this.genome = null;
    }

    public IGVBEDCodec(Genome genome) {
        this.genome = genome;
    }

    public void setGffTags(boolean gffTags) {
        this.gffTags = gffTags;
    }

    public boolean isGffTags() {
        return this.gffTags;
    }

    public BasicFeature decode(String nextLine) {

        if (nextLine.trim().length() == 0) {
            return null;
        }

        if (nextLine.startsWith("#") || nextLine.startsWith("track") || nextLine.startsWith("browser")) {
            this.readHeaderLine(nextLine);
            return null;
        }

        String [] tokens = Globals.whitespacePattern.split(nextLine);
        int tokenCount = tokens.length;

        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.

        if (tokenCount < 2) {
            return null;
        }

        String c = tokens[0];
        String chr = genome == null ? c : genome.getChromosomeAlias(c);

        //BED format, and IGV, use starting element as 0.
        int start = Integer.parseInt(tokens[1]) - startOffsetValue;

        int end = start + 1;
        if (tokenCount > 2) {
            end = Integer.parseInt(tokens[2]);
        }

        BasicFeature feature = spliceJunctions ?
                new SpliceJunctionFeature(chr, start, end) :
                new BasicFeature(chr, start, end);

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

                feature.setAttributes(atts);


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
                if (spliceJunctions) {
                    ((SpliceJunctionFeature) feature).setJunctionDepth((int) score);
                }
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

        // Color
        if (tokenCount > 8) {
            String colorString = tokens[8];
            if (colorString.trim().length() > 0 && !colorString.equals(".")) {
                feature.setColor(ParsingUtils.parseColor(colorString));
            }
        }

        // Exons
        if (tokenCount > 11) {
            createExons(start, tokens, feature, chr, feature.getStrand());
            //todo: some refactoring that allows this hack to be removed
            if (spliceJunctions) {
                SpliceJunctionFeature junctionFeature = (SpliceJunctionFeature) feature;

                List<Exon> exons = feature.getExons();

                junctionFeature.setJunctionStart(start + exons.get(0).getLength());
                junctionFeature.setJunctionEnd(end - exons.get(1).getLength());

            }
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

        int cdStart = Integer.parseInt(tokens[6]) - startOffsetValue;
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

    /**
     * Encode a feature as a BED string.
     *
     * @param feature - feature to encode
     * @return the encodec string
     */
    public String encode(BasicFeature feature) {

        StringBuffer buffer = new StringBuffer();

        buffer.append(feature.getChr());
        buffer.append("\t");
        final int featureStart = feature.getStart();
        buffer.append(String.valueOf(featureStart));
        buffer.append("\t");
        buffer.append(String.valueOf(feature.getEnd()));

        if (feature.getName() != null || (gffTags && feature.getDescription() != null)) {

            buffer.append("\t");

            if (gffTags && feature.getDescription() != null) {
                // mRNA<br>ID = LOC_Os01g01010.2<br>Name = LOC_Os01g01010.2<br>Parent = LOC_Os01g01010<br>
                //ID=LOC_Os01g01010.1:exon_1;Parent=LOC_Os01g01010.1
                String[] attrs = BR_PATTERN.split(feature.getDescription());
                buffer.append("\"");
                for (String att : attrs) {
                    String[] kv = EQ_PATTERN.split(att, 2);
                    if (kv.length > 1) {
                        buffer.append(kv[0].trim());
                        buffer.append("=");
                        String value = kv[1].trim();
                        buffer.append(URLEncoder.encode(value));
                        buffer.append(";");
                    }
                }
                buffer.append("\"");
            } else {
                buffer.append(feature.getName());
            }

            boolean more = !Float.isNaN(feature.getScore()) || feature.getStrand() != Strand.NONE ||
                    feature.getColor() != null || feature.getExonCount() > 0;

            if (more) {
                buffer.append("\t");
                // UCSC scores are integers between 0 and 1000, but
                float score = feature.getScore();
                if (Float.isNaN(score)) {
                    buffer.append("1000");

                } else {
                    boolean isInt = (Math.floor(score) == score);
                    buffer.append(String.valueOf(isInt ? (int) score : score));
                }


                more = feature.getStrand() != Strand.NONE || feature.getColor() != null || feature.getExonCount() > 0;
                if (more) {
                    buffer.append("\t");
                    Strand strand = feature.getStrand();
                    if (strand == Strand.NONE) buffer.append(" ");
                    else if (strand == Strand.POSITIVE) buffer.append("+");
                    else if (strand == Strand.NEGATIVE) buffer.append("-");

                    more = feature.getColor() != null || feature.getExonCount() > 0;

                    if (more) {
                        // Must continue if feature has color or exons
                        java.util.List<Exon> exons = feature.getExons();
                        if (feature.getColor() != null || exons != null) {
                            buffer.append("\t");
                            buffer.append(String.valueOf(feature.getThickStart()));
                            buffer.append("\t");
                            buffer.append(String.valueOf(feature.getThickEnd()));
                            buffer.append("\t");

                            java.awt.Color c = feature.getColor();
                            buffer.append(c == null ? "." : ColorUtilities.colorToString(c));
                            buffer.append("\t");

                            if (exons != null && exons.size() > 0) {
                                buffer.append(String.valueOf(exons.size()));
                                buffer.append("\t");

                                for (Exon exon : exons) {
                                    buffer.append(String.valueOf(exon.getLength()));
                                    buffer.append(",");
                                }
                                buffer.append("\t");
                                for (Exon exon : exons) {
                                    int exonStart = exon.getStart() - featureStart;
                                    buffer.append(String.valueOf(exonStart));
                                    buffer.append(",");
                                }

                            }
                        }
                    }
                }
            }
        }

        return buffer.toString();
    }


}


