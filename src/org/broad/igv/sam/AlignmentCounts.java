package org.broad.igv.sam;

import java.util.Set;

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

    int getMaxCount();

    String getValueStringAt(int pos);

    PositionIterator getPositionIterator();

    public boolean isMismatch(int pos, char ref, String chr, float snpThreshold);


    static interface PositionIterator {
        int nextPosition();
    }
}
