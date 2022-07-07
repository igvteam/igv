package org.broad.igv.sam;

import java.util.Arrays;

/**
 * A read-only implementation of a byte sub-array.  Created to support alignment blocks.
 */
public class ByteSubarray {

    public byte [] backingArray;
    public int startOffset;
    public int length;
    byte fillByte;

    public ByteSubarray(byte[] backingArray, int startOffset, int length, byte fillByte) {
        this.backingArray = backingArray;
        this.startOffset = startOffset;
        this.length = length;
        this.fillByte = fillByte;
    }

    public byte getByte(int idx) {

        if(idx < 0 || idx >= length) {
            throw new IndexOutOfBoundsException("Index out of bounds: "  + idx);
        }

        int i = startOffset + idx;
        return i < backingArray.length ?  backingArray[startOffset + idx] : fillByte;
    }

    public byte [] copyOfRange(int start, int end) {
        return Arrays.copyOfRange(backingArray, startOffset + start, startOffset + end);
    }

    public String getString() {
        return new String(backingArray, startOffset, length);
    }

    @Deprecated
    public byte [] getBytes() {
        return Arrays.copyOfRange(backingArray, startOffset, startOffset + length);
    }

}
