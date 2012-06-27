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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.dranger;

import org.apache.batik.dom.svg12.Global;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class DRangerParser {
    int numColumn;
    int chr1Column;
    int str1Column;
    int pos1Column;
    int chr2Column;
    int str2Column;
    int pos2Column;
    int tumreadsColumn;
    int normreadsColumn;
    int classColumn;
    int spanColumn;
    int site1Column;
    int site2Column;
    int qualityColumn;
    int scoreColumn;

    static Logger log = Logger.getLogger(DRangerParser.class);

    //num	chr1	str1	pos1	chr2	str2	pos2
    //1	    chr1	(+)	21083863	chr1	(+)	25279578
    // tumreads	normreads	class	span	site1	site2	quality	score
    // 48	    0	    long_range	4195715	3'-UTR of ENST0000037493
    public List<FeatureTrack> loadTracks(ResourceLocator locator, Genome genome) {

        List<FeatureTrack> tracks = new ArrayList();
        AsciiLineReader reader = null;
        List<org.broad.tribble.Feature> features = new ArrayList(5000);

        int parseColumn = -1;
        String nextLine = null;
        int rowCounter = 0;

        try {
            reader = ParsingUtils.openAsciiReader(locator);
            setColumns(reader.readLine());
            rowCounter++;
            while ((nextLine = reader.readLine()) != null) {

               String[] tokens = Globals.tabPattern.split(nextLine, -1);
               rowCounter++;

               int nTokens = tokens.length;
                if (nTokens > pos2Column) {
                    int index = Integer.parseInt(tokens[numColumn]);

                    String chr1 = genome.getChromosomeAlias(tokens[chr1Column]);
                    int pos1 = Integer.parseInt(tokens[pos1Column]);
                    String str1 = tokens[str1Column];
                    Strand strand1 = (str1.equals("0") || str1.equals("(+)")) ? Strand.POSITIVE : Strand.NEGATIVE;

                    String chr2 = genome.getChromosomeAlias(tokens[chr2Column]);
                    int pos2 = Integer.parseInt(tokens[pos2Column]);
                    String str2 = tokens[str2Column];
                    Strand strand2 = (str2.equals("0") || str2.equals("(+)")) ? Strand.POSITIVE : Strand.NEGATIVE;

                    DRangerFeature feature = new DRangerFeature(chr1, pos1, strand1, chr2, pos2, strand2);
                    feature.setIndex(index);

                    if (tumreadsColumn > 0) {
                        feature.setTumreads(toInt(tokens[tumreadsColumn]));
                    }
                    if (normreadsColumn > 0) {
                        feature.setNormreads(toInt(tokens[normreadsColumn]));
                    }
                    if (classColumn > 0) {
                        feature.setFeatureClass(tokens[classColumn]);
                    }
                    if (spanColumn > 0) {
                        feature.setSpan(toInt(tokens[spanColumn]));
                    }
                    if (site1Column > 0) {
                        feature.setSite1(tokens[site1Column]);
                    }
                    if (site2Column > 0) {
                        feature.setSite2(tokens[site2Column]);
                    }
                    if (qualityColumn > 0) {
                        feature.setQuality(toInt(tokens[qualityColumn]));
                    }
                    if (scoreColumn > 0) {
                        feature.setScore(toInt(tokens[scoreColumn]));
                    }

                    features.add(feature);
                }
            }
        } catch (NumberFormatException ne) {
            throw new ParserException("Column " + parseColumn + " must be a numeric value", rowCounter, nextLine);
        } catch (Exception e) {
            log.error("Error parsing dRanger file", e);
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }

        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features, genome));
        track.setName(locator.getTrackName());
        track.setRendererClass(DRangerRenderer.class);
        tracks.add(track);

        return tracks;
    }

    private int toInt(String token) {
        try {
            return Integer.parseInt(token);
        } catch (NumberFormatException e) {
            return 0;
        }
    }


    private void setColumns(String header) throws Exception {
        String[] tokens = header.split("\t");

        Map<String, Integer> map = new HashMap();
        for (int i = 0; i < tokens.length; i++) {
            map.put(tokens[i].toLowerCase(), i);
        }

        chr1Column = map.get("chr1").intValue();
        str1Column = map.get("str1").intValue();
        pos1Column = map.get("pos1").intValue();
        chr2Column = map.get("chr2").intValue();
        str2Column = map.get("str2").intValue();
        pos2Column = map.get("pos2").intValue();

        if (map.containsKey("num")) {
            numColumn = map.get("num").intValue();
        }
        if (map.containsKey("tumreads")) {
            tumreadsColumn = map.get("tumreads").intValue();
        }
        if (map.containsKey("normreads")) {
            tumreadsColumn = map.get("normreads").intValue();
        }
        if (map.containsKey("class")) {
            classColumn = map.get("class").intValue();
        }
        if (map.containsKey("span")) {
            spanColumn = map.get("span").intValue();
        }
        if (map.containsKey("site1")) {
            site1Column = map.get("site1").intValue();
        }
        if (map.containsKey("site2")) {
            site2Column = map.get("site2").intValue();
        }
        if (map.containsKey("quality")) {
            qualityColumn = map.get("quality").intValue();
        }
        if (map.containsKey("score")) {
            scoreColumn = map.get("score").intValue();
        }
    }
}
