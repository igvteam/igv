/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Fred Hutchinson Cancer Research Center and Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
                break;
            case REVERSE:
                junctions = new ArrayList<SpliceJunctionFeature>(negStartEndJunctionsMap.values());
                break;
            case BOTH:
                junctions = new ArrayList<SpliceJunctionFeature>(posStartEndJunctionsMap.values());
                junctions.addAll(negStartEndJunctionsMap.values());
                break;
            default:
                junctions = combineStrandJunctionsMaps();
        }

        List<SpliceJunctionFeature> filteredJunctions = filterJunctionList(this.loadOptions, junctions);

        FeatureUtils.sortFeatureList(filteredJunctions);

        return filteredJunctions;

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
        if (gapTypes != null) {
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


    /**
     * Combine junctions from both strands.  Used for Sashimi plot.
     * Note: Flanking depth arrays are not combined.
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
                // No existing (+) junction here, just add the (-) one\
                combinedStartEndJunctionsMap.put(junctionStart, junctionEnd, negFeat);
            } else {
                int newJunctionDepth = junction.getJunctionDepth() + negFeat.getJunctionDepth();
                junction.setJunctionDepth(newJunctionDepth);
            }
        }

        return new ArrayList<SpliceJunctionFeature>(combinedStartEndJunctionsMap.values());
    }


    void setLoadOptions(LoadOptions loadOptions) {
        int oldMinJunctionCoverage = this.loadOptions.minJunctionCoverage;
        //Can't change this, need to reload everything
        assert this.loadOptions.minReadFlankingWidth == loadOptions.minReadFlankingWidth;
        this.loadOptions = loadOptions;

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
