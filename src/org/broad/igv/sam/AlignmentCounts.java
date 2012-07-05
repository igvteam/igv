package org.broad.igv.sam;

/**
 * @author Jim Robinson
 * @date 11/22/11
 */
public interface AlignmentCounts {

    int getTotalCount(int pos);

    int getNegTotal(int pos);

    int getPosTotal(int pos);

    int getTotalQuality(int pos);

    int getCount(int pos, byte b);

    int getNegCount(int pos, byte b);

    int getPosCount(int pos, byte b);

    int getDelCount(int pos);

    int getInsCount(int pos);

    int getQuality(int pos, byte b);

    int getAvgQuality(int pos, byte b);

    void incCounts(Alignment alignment);

    int getStart();

    int getEnd();

    int getNumberOfPoints();

    int getMaxCount();

    String getValueStringAt(int pos);

    public boolean isMismatch(int pos, byte ref, String chr, float snpThreshold);

    BisulfiteCounts getBisulfiteCounts();

    void finish();


    static interface PositionIterator {
        int nextPosition();
    }
}
