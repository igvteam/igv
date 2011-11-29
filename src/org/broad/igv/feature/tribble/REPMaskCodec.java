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

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;

/**
 * Basically BED format with some columns rearranged
 */
public class REPMaskCodec implements FeatureCodec {

    // Declare a  array once, to be reused.
    String[] tokens = new String[15];
    int startBase = 0;
    FeatureFileHeader header;

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
        }
        catch (IOException e) {
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
        return true;
    }

    public Feature decodeLoc(String line) {
        return decode(line);
    }

    public BasicFeature decode(String nextLine) {

        if (nextLine.trim().length() == 0 || nextLine.startsWith("#")) {
            return null;
        }

        int tokenCount = ParsingUtils.splitWhitespace(nextLine, tokens);

        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.

        if (tokenCount < 2) {
            return null;
        }

        String chr = tokens[0];  //genome == null ? tokens[0] : genome.getChromosomeAlias(tokens[0]);
        int start = Integer.parseInt(tokens[1]) - startBase;

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
            char strand = (strandString.length() == 0)  ? ' ' : strandString.charAt(0);

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

    public Class getFeatureType() {
        return BasicFeature.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

}