/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */


package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.Strand;

import java.io.IOException;
import java.util.*;

/**
 * A helper class for computing splice junctions from alignments.
 * <p/>
 * dhmay 20111014 moving min junction coverage and min alignment flanking width references to preferences
 *
 * @author dhmay, jrobinso
 * @date Jul 3, 2011
 */
public class SpliceJunctionHelper {

    static Logger log = Logger.getLogger(SpliceJunctionHelper.class);


    public static List<SpliceJunctionFeature> computeFeatures(Iterator<Alignment> iterator) throws IOException {

        List<SpliceJunctionFeature> spliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>();


        //we need to keep the positive and negative strand junctions separate, since
        //they don't represent the same thing and are rendered separately
        Map<Integer, Map<Integer, SpliceJunctionFeature>> posStartEndJunctionsMap = new HashMap<Integer, Map<Integer, SpliceJunctionFeature>>();
        Map<Integer, Map<Integer, SpliceJunctionFeature>> negStartEndJunctionsMap = new HashMap<Integer, Map<Integer, SpliceJunctionFeature>>();

        //dhmay adding after this was moved from track level to global
        PreferenceManager prefs = PreferenceManager.getInstance();
        int minJunctionCoverage = prefs.getAsInt(PreferenceManager.SAM_JUNCTION_MIN_COVERAGE);
        int minReadFlankingWidth = prefs.getAsInt(PreferenceManager.SAM_JUNCTION_MIN_FLANKING_WIDTH);


        while (iterator.hasNext()) {
            //Any alignment with 2 or more blocks is considered to be a splice junction
            Alignment alignmentFromIterator = iterator.next();

            //This is an odd bit of code.  If alignmentFromIterator is a paired alignment, then we actually need
            //to do the processing separately on each of the paired alignments.  Otherwise, we need to process
            //alignmentFromIterator, itself.
            List<Alignment> alignmentsRepresentedByThisAlignment = new ArrayList<Alignment>();
            if (alignmentFromIterator instanceof PairedAlignment) {
                PairedAlignment alAsPair = (PairedAlignment) alignmentFromIterator;
                if (alAsPair.getFirstAlignment() != null) {
                    alignmentsRepresentedByThisAlignment.add(alAsPair.getFirstAlignment());
                }
                if (alAsPair.getSecondAlignment() != null) {
                    alignmentsRepresentedByThisAlignment.add(alAsPair.getSecondAlignment());
                }
            } else alignmentsRepresentedByThisAlignment.add(alignmentFromIterator);

            for (Alignment alignment : alignmentsRepresentedByThisAlignment) {
                AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                Object strandAttr = alignment.getAttribute("XS");
                if (blocks.length < 2 || strandAttr == null) {
                    continue;
                }

                //there may be other ways in which this is indicated. May have to code for them later
                boolean isNegativeStrand = strandAttr.toString().charAt(0) == '-';

                Map<Integer, Map<Integer, SpliceJunctionFeature>> startEndJunctionsMapThisStrand =
                        isNegativeStrand ? negStartEndJunctionsMap : posStartEndJunctionsMap;

                int flankingStart = -1;
                int junctionStart = -1;
                int gapCount = -1;
                char[] gapTypes = alignment.getGapTypes();
                //for each pair of blocks, create or add evidence to a splice junction
                for (AlignmentBlock block : blocks) {
                    int flankingEnd = block.getEnd();
                    int junctionEnd = block.getStart();
                    if (junctionStart != -1 && gapCount < gapTypes.length && gapTypes[gapCount] == SamAlignment.SKIPPED_REGION) {
                        //only proceed if the flanking regions are both bigger than the minimum
                        if (minReadFlankingWidth == 0 ||
                                ((junctionStart - flankingStart >= minReadFlankingWidth) &&
                                        (flankingEnd - junctionEnd >= minReadFlankingWidth))) {
                            Map<Integer, SpliceJunctionFeature> endJunctionsMap =
                                    startEndJunctionsMapThisStrand.get(junctionStart);
                            if (endJunctionsMap == null) {
                                endJunctionsMap = new HashMap<Integer, SpliceJunctionFeature>();
                                startEndJunctionsMapThisStrand.put(junctionStart, endJunctionsMap);
                            }
                            SpliceJunctionFeature junction = endJunctionsMap.get(junctionEnd);
                            if (junction == null) {
                                junction = new SpliceJunctionFeature(alignment.getChr(), junctionStart, junctionEnd,
                                        isNegativeStrand ? Strand.NEGATIVE : Strand.POSITIVE);
                                endJunctionsMap.put(junctionEnd, junction);
                                spliceJunctionFeatures.add(junction);
                            }
                            junction.addRead(flankingStart, flankingEnd);
                        }
                    }
                    flankingStart = junctionEnd;
                    junctionStart = flankingEnd;
                    gapCount += 1;
                }
            }
        }

        //get rid of any features without enough coverage
        if (minJunctionCoverage > 1) {
            List<SpliceJunctionFeature> coveredFeatures =
                    new ArrayList<SpliceJunctionFeature>(spliceJunctionFeatures.size());
            for (SpliceJunctionFeature feature : spliceJunctionFeatures)
                if (feature.getJunctionDepth() >= minJunctionCoverage)
                    coveredFeatures.add(feature);
            spliceJunctionFeatures = coveredFeatures;
        }

        //Sort by increasing beginning of start flanking region, as required by the renderer
        Collections.sort(spliceJunctionFeatures, new Comparator<IGVFeature>() {
            public int compare(IGVFeature o1, IGVFeature o2) {
                return o1.getStart() - o2.getStart();
            }
        });


        return spliceJunctionFeatures;

    }

}
