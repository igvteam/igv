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

package org.broad.igv.feature.tribble;

import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.collections.MultiMap;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Basically BED format with some columns rearranged
 * <p/> Columns, from UCSC documentation
 * <p/>
 * 0  bin	585	smallint(5) unsigned	Indexing field to speed chromosome range queries.
 * 1  swScore	1504	int(10) unsigned	Smith Waterman alignment score
 * 2  milliDiv	13	int(10) unsigned	Base mismatches in parts per thousand
 * 3  milliDel	4	int(10) unsigned	Bases deleted in parts per thousand
 * 4  milliIns	13	int(10) unsigned	Bases inserted in parts per thousand
 * 5  genoName	chr1	varchar(255)	Genomic sequence name
 * 6  genoStart	10000	int(10) unsigned	Start in genomic sequence
 * 7  genoEnd	10468	int(10) unsigned	End in genomic sequence
 * 8  genoLeft	-249240153	int(11)	-#bases after match in genomic sequence
 * 9  strand	+	char(1)	Relative orientation + or -
 * 10 repName	(CCCTAA)n	varchar(255)	Name of repeat
 * 11 repClass	Simple_repeat	varchar(255)	Class of repeat
 * 12 repFamily	Simple_repeat	varchar(255)	Family of repeat
 * 13 repStart	1	int(11)	Start (if strand is +) or -#bases after match (if strand is -) in repeat sequence
 * 14 repEnd	463	int(11)	End in repeat sequence
 * 15 repLeft	0	int(11)	-#bases after match (if strand is +) or start (if strand is -) in repeat sequence
 * 16 id	1	char(1)	First digit of id field in RepeatMasker .out file. Best ignored.
 */
public class REPMaskCodec extends AsciiFeatureCodec<BasicFeature> {

    FeatureFileHeader header;
    Genome genome;

    public REPMaskCodec(Genome genome) {
        super(BasicFeature.class);
        this.genome = genome;
    }

    public Object readHeader(LineReader reader) {

        String nextLine;
        header = new FeatureFileHeader();
        header.setTrackType(TrackType.REPMASK);
        int nLines = 0;

        try {
            while ((nextLine = reader.readLine()) != null &&
                    (nextLine.startsWith("#") || nextLine.startsWith("track")) ||
                    nextLine.startsWith("browser")) {
                nLines++;
                if (nextLine.startsWith("#type")) {
                    String[] tokens = nextLine.split("=");
                    if (tokens.length > 1) {
                        try {
                            header.setTrackType(TrackType.valueOf(tokens[1]));
                        } catch (Exception e) {
                            // log.error("Error converting track type: " + tokens[1]);
                        }
                    }
                } else if (nextLine.startsWith("track")) {
                    TrackProperties tp = new TrackProperties();
                    ParsingUtils.parseTrackLine(nextLine, tp);
                    header.setTrackProperties(tp);
                }
            }
            return header;
        } catch (IOException e) {
            throw new CodecLineParsingException("Error parsing header", e);
        }
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
        return true; // Optimisitic!
    }

    /**
     * Return an abbreviated feature, containing only location information.  Used for indexing.
     *
     * @param line
     * @return
     */
    public Feature decodeLoc(String line) {
        String[] tokens = Globals.singleTabMultiSpacePattern.split(line);
        if(tokens.length < 15) {
            return decodeLegacy(tokens);
        }
        String chr = genome == null ? tokens[5] : genome.getChromosomeAlias(tokens[5]);
        int start = Integer.parseInt(tokens[6]);
        int end = Integer.parseInt(tokens[7]);
        return new BasicFeature(chr, start, end);
    }

    public BasicFeature decode(String nextLine) {

        if (nextLine.trim().length() == 0 || nextLine.startsWith("#")) {
            return null;
        }

        String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
        int tokenCount = tokens.length;

        if (tokenCount < 15) {
            return decodeLegacy(tokens);
        }

        String chr = genome == null ? tokens[5] : genome.getChromosomeAlias(tokens[5]);
        int start = Integer.parseInt(tokens[6]);
        int end = Integer.parseInt(tokens[7]);
        BasicFeature feature = new BasicFeature(chr, start, end);

        String strandString = tokens[3].trim();
        char strand = (strandString.length() == 0) ? ' ' : strandString.charAt(0);
        if (strand == '-') {
            feature.setStrand(Strand.NEGATIVE);
        } else if (strand == '+') {
            feature.setStrand(Strand.POSITIVE);
        } else {
            feature.setStrand(Strand.NONE);
        }
        String name = tokens[10];
        feature.setName(name);
        feature.setIdentifier(name);

        MultiMap<String, String> attributes = new MultiMap<String, String>();
        attributes.put("Smith Waterman score", tokens[1]);
        attributes.put("base mismatches per thousand", tokens[2]);
        attributes.put("bases deleted per thousand", tokens[3]);
        attributes.put("bases inserted per thousand", tokens[4]);
        attributes.put("repeat class", tokens[11]);
        attributes.put("repeat family", tokens[12]);
        attributes.put("repeat start", tokens[13]);
        attributes.put("repeat end", tokens[14]);
        feature.setAttributes(attributes);

        return feature;
    }


    public BasicFeature decodeLegacy(String [] tokens) {


        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.
        int tokenCount = tokens.length;

        if (tokenCount < 2) {
            return null;
        }

        String chr = genome == null ? tokens[0] : genome.getChromosomeAlias(tokens[0]);
        int start = Integer.parseInt(tokens[1]);

        int end = start + 1;
        if (tokenCount > 2) {
            end = Integer.parseInt(tokens[2]);
        }

        BasicFeature feature = new BasicFeature(chr, start, end);

        // The rest of the columns are optional.  Stop parsing upon encountering
        // a non-expected value

        // Strand
        if (tokenCount > 3) {
            String strandString = tokens[3].trim();
            char strand = (strandString.length() == 0) ? ' ' : strandString.charAt(0);

            if (strand == '-') {
                feature.setStrand(Strand.NEGATIVE);
            } else if (strand == '+') {
                feature.setStrand(Strand.POSITIVE);
            } else {
                feature.setStrand(Strand.NONE);
            }
        }

        // Name
        if (tokenCount > 4) {
            String name = tokens[4].replaceAll("\"", "");
            feature.setName(name);
            feature.setIdentifier(name);
        }


        return feature;
    }

}
