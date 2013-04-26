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

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.Strand;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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

    List<SpliceJunctionFeature> spliceJunctionFeatures = new ArrayList();

    Table<Integer, Integer, SpliceJunctionFeature> posStartEndJunctionsMap = HashBasedTable.create();
    Table<Integer, Integer, SpliceJunctionFeature> negStartEndJunctionsMap = HashBasedTable.create();

    private final LoadOptions loadOptions;

    public SpliceJunctionHelper(LoadOptions loadOptions){
        this.loadOptions = loadOptions;
    }

    public List<SpliceJunctionFeature> getFeatures() {
        return spliceJunctionFeatures;

    }

    public void addAlignment(Alignment alignment) {

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        if (blocks == null || blocks.length < 2) {
            return;
        }

        //there may be other ways in which this is indicated. May have to code for them later
        boolean isNegativeStrand;
        Object strandAttr = alignment.getAttribute("XS");
        if (strandAttr != null) {
            isNegativeStrand = strandAttr.toString().charAt(0) == '-';
        } else {
            isNegativeStrand = alignment.isNegativeStrand(); // <= TODO -- this isn't correct for all libraries.
        }

        Table<Integer, Integer, SpliceJunctionFeature> startEndJunctionsTableThisStrand =
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
                if (loadOptions.minReadFlankingWidth == 0 ||
                        ((junctionStart - flankingStart >= loadOptions.minReadFlankingWidth) &&
                                (flankingEnd - junctionEnd >= loadOptions.minReadFlankingWidth))) {

                    SpliceJunctionFeature junction = startEndJunctionsTableThisStrand.get(junctionStart, junctionEnd);

                    if (junction == null) {
                        junction = new SpliceJunctionFeature(alignment.getChr(), junctionStart, junctionEnd,
                                isNegativeStrand ? Strand.NEGATIVE : Strand.POSITIVE);
                        startEndJunctionsTableThisStrand.put(junctionStart, junctionEnd, junction);
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


    public void finish() {
        //get rid of any features without enough coverage
        if (loadOptions.minJunctionCoverage > 1) {
            List<SpliceJunctionFeature> coveredFeatures = new ArrayList<SpliceJunctionFeature>(spliceJunctionFeatures.size());
            for (SpliceJunctionFeature feature : spliceJunctionFeatures) {
                if (feature.getJunctionDepth() >= loadOptions.minJunctionCoverage) {
                    coveredFeatures.add(feature);
                }
            }
            spliceJunctionFeatures = coveredFeatures;
        }

        //Sort by increasing beginning of start flanking region, as required by the renderer
        Collections.sort(spliceJunctionFeatures, new Comparator<IGVFeature>() {
            public int compare(IGVFeature o1, IGVFeature o2) {
                return o1.getStart() - o2.getStart();
            }
        });
    }

    public static class LoadOptions {

        private static PreferenceManager prefs = PreferenceManager.getInstance();

        public final boolean showSpliceJunctions;
        public final int minJunctionCoverage;
        public final int minReadFlankingWidth;

        public LoadOptions(boolean showSpliceJunctions){
            this(showSpliceJunctions, prefs.getAsInt(PreferenceManager.SAM_JUNCTION_MIN_COVERAGE),
                    prefs.getAsInt(PreferenceManager.SAM_JUNCTION_MIN_FLANKING_WIDTH));
        }

        public LoadOptions(boolean showSpliceJunctions, int minJunctionCoverage, int minReadFlankingWidth){
            this.showSpliceJunctions = showSpliceJunctions;
            this.minJunctionCoverage = minJunctionCoverage;
            this.minReadFlankingWidth = minReadFlankingWidth;
        }
    }

}
