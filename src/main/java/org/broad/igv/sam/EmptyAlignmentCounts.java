package org.broad.igv.sam;

public class EmptyAlignmentCounts implements AlignmentCounts {

    private static EmptyAlignmentCounts instance;

    public static AlignmentCounts getInstance() {
        if(instance == null) {
            instance = new EmptyAlignmentCounts();
        }
        return instance;
    }

    @Override
    public void incCounts(Alignment alignment) {

    }

    @Override
    public int getTotalCount(int pos) {
        return 0;
    }

    @Override
    public int getTotalPositiveCount(int pos) {
        return 0;
    }

    @Override
    public int getTotalNegativeCount(int pos) {
        return 0;
    }

    @Override
    public int getTotalQuality(int pos) {
        return 0;
    }

    @Override
    public int getCount(int pos, byte b) {
        return 0;
    }

    @Override
    public int getNegCount(int pos, byte b) {
        return 0;
    }

    @Override
    public int getPosCount(int pos, byte b) {
        return 0;
    }

    @Override
    public int getDelCount(int pos) {
        return 0;
    }

    @Override
    public int getInsCount(int pos) {
        return 0;
    }

    @Override
    public int getQuality(int pos, byte b) {
        return 0;
    }

    @Override
    public int getNumberOfPoints() {
        return 0;
    }

    @Override
    public int getMaxCount(int origin, int end) {
        return 0;
    }

    @Override
    public String getValueStringAt(int pos) {
        return null;
    }

    @Override
    public boolean isConsensusMismatch(int pos, byte ref, String chr, float snpThreshold) {
        return false;
    }

    @Override
    public boolean isConsensusDeletion(int start, int end, float snpThreshold) {
        return false;
    }

    @Override
    public boolean isConsensusInsertion(int pos, float snpThreshold) {
        return false;
    }

    @Override
    public BisulfiteCounts getBisulfiteCounts() {
        return null;
    }

    @Override
    public int getBucketSize() {
        return 0;
    }

    @Override
    public boolean hasBaseCounts() {
        return false;
    }

    @Override
    public void finish() {

    }

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }
}
