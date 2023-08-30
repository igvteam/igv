package org.broad.igv.feature.tribble.reader;


public class VPointer {

    private static final int SHIFT_AMOUNT = 16;
    private static final int OFFSET_MASK = 0xffff;
    private static final long ADDRESS_MASK = 0xFFFFFFFFFFFFL;


    long block;
    int offset;

    public VPointer(long v) {
        this.block = (v >> SHIFT_AMOUNT) & ADDRESS_MASK;
        this.offset = (int) (v & OFFSET_MASK);
    }

    public boolean isLessThan(VPointer vp) {
        return this.block < vp.block ||
                (this.block == vp.block && this.offset < vp.offset);
    }

    public boolean isGreaterThan(VPointer vp) {
        return this.block > vp.block ||
                (this.block == vp.block && this.offset > vp.offset);
    }

}