package org.igv.sam;

import htsjdk.tribble.Feature;
import org.igv.sam.mods.BaseModificationCounts;

/**
 * @author Jim Robinson
 * @date 11/22/11
 */
public interface AlignmentCounts extends Feature {

    void incCounts(Alignment alignment);

    int getTotalCount(int pos);

    public int getTotalPositiveCount(int pos);

    public int getTotalNegativeCount(int pos);

    int getTotalQuality(int pos);

    int getCount(int pos, byte b);

    int getNegCount(int pos, byte b);

    int getPosCount(int pos, byte b);

    int getDelCount(int pos);

    int getInsCount(int pos);

    int getQuality(int pos, byte b);

    int getNumberOfPoints();

    int getMaxCount(int origin, int end);

    String getValueStringAt(String chr, int pos);

    boolean isConsensusMismatch(int pos, byte ref, String chr, float snpThreshold);

    boolean isConsensusDeletion(int start, int end, float snpThreshold);

    boolean isConsensusInsertion(int pos, float snpThreshold);

    BisulfiteCounts getBisulfiteCounts();

    default BaseModificationCounts getModifiedBaseCounts() {
        return null;
    }

    int getBucketSize();

    boolean hasBaseCounts();

    void finish();

}
