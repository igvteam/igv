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
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.Strand;

import java.util.ArrayList;
import java.util.List;

/**
 * A helper class for computing splice junctions from alignments.
 * Junctions are filtered based on minimum flanking width on loading, so data
 * needs to be
 *
 * @author dhmay, jrobinso
 * @date Jul 3, 2011
 */
public class SpliceJunctionHelper {

    static Logger log = Logger.getLogger(SpliceJunctionHelper.class);

    List<SpliceJunctionFeature> allSpliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>();
    //  List<SpliceJunctionFeature> filteredSpliceJunctionFeatures = null;
    List<SpliceJunctionFeature> filteredCombinedFeatures = null;

    Table<Integer, Integer, SpliceJunctionFeature> posStartEndJunctionsMap = HashBasedTable.create();
    Table<Integer, Integer, SpliceJunctionFeature> negStartEndJunctionsMap = HashBasedTable.create();

    private LoadOptions loadOptions;

    public SpliceJunctionHelper(LoadOptions loadOptions) {
        this.loadOptions = loadOptions;
    }

    public List<SpliceJunctionFeature> getFilteredJunctions(SpliceJunctionFinderTrack.StrandOption strandOption) {

        List<SpliceJunctionFeature> junctions;

        switch (strandOption) {
            case FORWARD:
                junctions = new ArrayList<SpliceJunctionFeature>(posStartEndJunctionsMap.values());
            case REVERSE:
                junctions = new ArrayList<SpliceJunctionFeature>(negStartEndJunctionsMap.values());
            case IGNORE:
                junctions = combineStrandJunctionsMaps();
            default:
                junctions = new ArrayList<SpliceJunctionFeature>(posStartEndJunctionsMap.values());
                junctions.addAll(negStartEndJunctionsMap.values());
        }

        filterJunctionList(this.loadOptions, junctions);

        FeatureUtils.sortFeatureList(junctions);

        return junctions;

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
            if (alignment.isPaired()) {
                isNegativeStrand = alignment.getFirstOfPairStrand() == Strand.NEGATIVE;
            } else {
                isNegativeStrand = alignment.isNegativeStrand(); // <= TODO -- this isn't correct for all libraries.
            }
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
            if (junctionStart != -1 && gapCount < gapTypes.length && gapTypes[gapCount] == SAMAlignment.SKIPPED_REGION) {

                //only proceed if the flanking regions are both bigger than the minimum
                if (loadOptions.minReadFlankingWidth == 0 ||
                        ((junctionStart - flankingStart >= loadOptions.minReadFlankingWidth) &&
                                (flankingEnd - junctionEnd >= loadOptions.minReadFlankingWidth))) {

                    SpliceJunctionFeature junction = startEndJunctionsTableThisStrand.get(junctionStart, junctionEnd);
                    if (junction == null) {
                        junction = new SpliceJunctionFeature(alignment.getChr(), junctionStart, junctionEnd,
                                isNegativeStrand ? Strand.NEGATIVE : Strand.POSITIVE);
                        startEndJunctionsTableThisStrand.put(junctionStart, junctionEnd, junction);
                        allSpliceJunctionFeatures.add(junction);
                    }
                    junction.addRead(flankingStart, flankingEnd);
                }

            }
            flankingStart = junctionEnd;
            junctionStart = flankingEnd;
            gapCount += 1;
        }
    }

    private static List<SpliceJunctionFeature> filterJunctionList(LoadOptions loadOptions, List<SpliceJunctionFeature> unfiltered) {

        if (loadOptions.minJunctionCoverage > 1) {
            List<SpliceJunctionFeature> coveredFeatures = new ArrayList<SpliceJunctionFeature>(unfiltered.size());
            for (SpliceJunctionFeature feature : unfiltered) {
                if (feature.getJunctionDepth() >= loadOptions.minJunctionCoverage) {
                    coveredFeatures.add(feature);
                }
            }
            return coveredFeatures;
        } else {
            return unfiltered;
        }
    }


    void setLoadOptions(LoadOptions loadOptions) {
        int oldMinJunctionCoverage = this.loadOptions.minJunctionCoverage;
        //Can't change this, need to reload everything
        assert this.loadOptions.minReadFlankingWidth == loadOptions.minReadFlankingWidth;
        this.loadOptions = loadOptions;

    }


    /**
     * We keep separate splice junction information by strand.
     * This combines both strand information.  Used for Sashimi plot,
     */
    private List<SpliceJunctionFeature> combineStrandJunctionsMaps() {

        // Start with all + junctions
        Table<Integer, Integer, SpliceJunctionFeature> combinedStartEndJunctionsMap = HashBasedTable.create(posStartEndJunctionsMap);

        // Merge in - junctions
        for (Table.Cell<Integer, Integer, SpliceJunctionFeature> negJunctionCell : negStartEndJunctionsMap.cellSet()) {

            int junctionStart = negJunctionCell.getRowKey();
            int junctionEnd = negJunctionCell.getColumnKey();
            SpliceJunctionFeature negFeat = negJunctionCell.getValue();

            SpliceJunctionFeature junction = combinedStartEndJunctionsMap.get(junctionStart, junctionEnd);

            if (junction == null) {
                // No existing (+) junction here, just add the (-) one
                junction = new SpliceJunctionFeature(negFeat.getChr(), junctionStart, junctionEnd, Strand.NONE);
                combinedStartEndJunctionsMap.put(junctionStart, junctionEnd, junction);
            } else {
                int newJunctionDepth = junction.getJunctionDepth() + negFeat.getJunctionDepth();
                junction.setJunctionDepth(newJunctionDepth);
            }
        }

        return new ArrayList<SpliceJunctionFeature>(combinedStartEndJunctionsMap.values());
    }

    public static class LoadOptions {

        private static PreferenceManager prefs = PreferenceManager.getInstance();

        public final int minJunctionCoverage;
        public final int minReadFlankingWidth;

        public LoadOptions() {
            this(prefs.getAsInt(PreferenceManager.SAM_JUNCTION_MIN_COVERAGE),
                    prefs.getAsInt(PreferenceManager.SAM_JUNCTION_MIN_FLANKING_WIDTH));
        }

        public LoadOptions(int minJunctionCoverage, int minReadFlankingWidth) {
            this.minJunctionCoverage = minJunctionCoverage;
            this.minReadFlankingWidth = minReadFlankingWidth;
        }
    }

}
