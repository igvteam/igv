package org.igv.sam;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import org.igv.logging.*;
import org.igv.feature.FeatureUtils;
import org.igv.feature.SpliceJunctionFeature;
import org.igv.feature.Strand;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;

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

    static Logger log = LogManager.getLogger(SpliceJunctionHelper.class);

    List<SpliceJunctionFeature> allSpliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>();
    //  List<SpliceJunctionFeature> filteredSpliceJunctionFeatures = null;
    List<SpliceJunctionFeature> filteredCombinedFeatures = null;

    Table<Integer, Integer, SpliceJunctionFeature> posStartEndJunctionsMap = HashBasedTable.create();
    Table<Integer, Integer, SpliceJunctionFeature> negStartEndJunctionsMap = HashBasedTable.create();


    public SpliceJunctionHelper() {
    }

    public List<SpliceJunctionFeature> getFilteredJunctions(SpliceJunctionTrack.StrandOption strandOption, int minJunctionCoverage) {

        List<SpliceJunctionFeature> junctions;

        switch (strandOption) {
            case FORWARD:
                junctions = new ArrayList<>(posStartEndJunctionsMap.values());
                break;
            case REVERSE:
                junctions = new ArrayList<>(negStartEndJunctionsMap.values());
                break;
            case BOTH:
                junctions = new ArrayList<>(posStartEndJunctionsMap.values());
                junctions.addAll(negStartEndJunctionsMap.values());
                break;
            default:
                junctions = combineStrandJunctionsMaps();
        }

        List<SpliceJunctionFeature> unfiltered = junctions;
        List<SpliceJunctionFeature> filteredJunctions;

        if (minJunctionCoverage > 1) {
            List<SpliceJunctionFeature> coveredFeatures = new ArrayList<SpliceJunctionFeature>(unfiltered.size());
            for (SpliceJunctionFeature feature : unfiltered) {
                if (feature.getJunctionDepth() >= minJunctionCoverage) {
                    coveredFeatures.add(feature);
                }
            }
            filteredJunctions = coveredFeatures;
        } else {
            filteredJunctions = unfiltered;
        }

        FeatureUtils.sortFeatureList(filteredJunctions);

        return filteredJunctions;

    }

    public void addAlignment(Alignment alignment) {

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        if (blocks == null || blocks.length < 2) {
            return;
        }

        // Determine strand.  First check for explicit attribute.
        boolean isNegativeStrand;
        Object strandAttr = alignment.getAttribute("TS");
        if(strandAttr == null) {
            strandAttr = alignment.getAttribute("XS");   // Older convention
        }

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


        // For each gap marked "skip" (cigar N), create or add evidence to a splice junction
        List<Gap> gaps = alignment.getGaps();
        int minReadFlankingWidth = PreferencesManager.getPreferences(Constants.RNA).getAsInt(Constants.SAM_JUNCTION_MIN_FLANKING_WIDTH);
        if (gaps != null) {
            for (Gap gap : gaps) {

                if (gap instanceof SpliceGap) {

                    SpliceGap spliceGap = (SpliceGap) gap;
                    //only proceed if the flanking regions are both bigger than the minimum
                    if (minReadFlankingWidth == 0 ||
                            (spliceGap.getFlankingLeft() >= minReadFlankingWidth &&
                                    spliceGap.getFlankingRight() >= minReadFlankingWidth)) {

                        int junctionStart = spliceGap.getStart();
                        int junctionEnd = junctionStart + spliceGap.getnBases();
                        int flankingStart = junctionStart - spliceGap.getFlankingLeft();
                        int flankingEnd = junctionEnd + spliceGap.getFlankingRight();

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
            }
        }
    }


    /**
     * Combine junctions from both strands.  Used for Sashimi plot.
     * Note: Flanking depth arrays are not combined.
     */
    private List<SpliceJunctionFeature> combineStrandJunctionsMaps() {

        // Start with all + junctions
        Table<Integer, Integer, SpliceJunctionFeature> combinedMap = HashBasedTable.create();

        for (SpliceJunctionFeature posFeature : posStartEndJunctionsMap.values()) {
            final int junctionStart = posFeature.getJunctionStart();
            final int junctionEnd = posFeature.getJunctionEnd();

            SpliceJunctionFeature combinedFeature = new SpliceJunctionFeature(posFeature.getChr(), junctionStart, junctionEnd);
            combinedFeature.setJunctionDepth(posFeature.getJunctionDepth());
            combinedMap.put(junctionStart, junctionEnd, combinedFeature);
        }


        // Merge in - junctions
        for (SpliceJunctionFeature negFeature : negStartEndJunctionsMap.values()) {

            int junctionStart = negFeature.getJunctionStart();
            int junctionEnd = negFeature.getJunctionEnd();

            SpliceJunctionFeature junction = combinedMap.get(junctionStart, junctionEnd);

            if (junction == null) {
                // No existing (+) junction here, just add the (-) one\
                SpliceJunctionFeature combinedFeature = new SpliceJunctionFeature(negFeature.getChr(), junctionStart, junctionEnd);
                combinedFeature.setJunctionDepth(negFeature.getJunctionDepth());
                combinedMap.put(junctionStart, junctionEnd, negFeature);
            } else {
                int newJunctionDepth = junction.getJunctionDepth() + negFeature.getJunctionDepth();
                junction.setJunctionDepth(newJunctionDepth);
            }
        }

        return new ArrayList<>(combinedMap.values());
    }


}
