package org.broad.igv.sam;

import java.util.Arrays;

/**
 * A read-only implementation of a byte sub-array.  Created to support aligment blocks.
 */
public class ByteSubarray {

    public byte [] backingArray;
    public int startOffset;
    public int length;

    public ByteSubarray(byte[] backingArray, int startOffset, int length) {
        this.backingArray = backingArray;
        this.startOffset = startOffset;
        this.length = length;
    }

    public byte getByte(int idx) {

        if(idx < 0 || idx >= length) {
            throw new IndexOutOfBoundsException("Index out of bounds: "  + idx);
        }

        return backingArray[startOffset + idx];
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
