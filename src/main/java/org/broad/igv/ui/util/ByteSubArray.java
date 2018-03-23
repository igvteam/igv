/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 UC San Diego
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

package org.broad.igv.ui.util;

/**
 * Created by jrobinson on 10/24/15.
 * <p/>
 * Class used to create a sub-array without making a copy of the original.  Initial use is for BAM alignment records
 */

public class ByteSubArray {

    byte[] bytes;
    short offset;
    public short length;

    public ByteSubArray(byte[] bytes, short offset, short length) {

        if (offset < 0 || offset + length > bytes.length) {
            throw new IndexOutOfBoundsException("Index out of bounds");
        }

        this.bytes = bytes;
        this.offset = offset;
        this.length = length;
    }

    public byte get(int i) {

        if (i < offset || i >= offset + length) {
            throw new IndexOutOfBoundsException();
        }

        return bytes[offset + i];
    }

    /**
     * Method provided to ease porting of old code.   This should not be used often as it defeats the purpose of
     * this class.
     *
     * @return
     */
    public byte[] getBytes() {
        byte[] subarray = new byte[length];
        System.arraycopy(bytes, offset, subarray, 0, length);
        return subarray;

    }
}
