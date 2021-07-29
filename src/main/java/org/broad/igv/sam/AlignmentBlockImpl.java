/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

public class AlignmentBlockImpl implements AlignmentBlock {

    private static ByteSubarray EMPTY_ARRAY = new ByteSubarray(new byte[0], 0, 0, (byte) 0);
    private int start;
    private int offset;
    private ByteSubarray bases;
    private int basesLength = -1;
    public ByteSubarray qualities;
    private boolean softClipped = false;
    private int pixelStart;
    private int pixelEnd;
    private int padding = 0;
    private char cigarOperator;

    public AlignmentBlockImpl(int start, byte[] bases, byte[] qualities) {
        this(start, bases, qualities, 0, bases.length, (char) 0);
    }

    public AlignmentBlockImpl(int start, byte[] bases, byte[] qualities, int offset, int nBases, char cigarOperator) {

        this.start = start;
        this.offset = offset;
        this.bases = new ByteSubarray(bases, offset, nBases, (byte) '?');
        this.basesLength = nBases;

        // qualities are optional in a SAMRecord, we might get null or an array of zero
        this.qualities = qualities == null || qualities.length == 0 ? EMPTY_ARRAY :  new ByteSubarray(qualities, offset, nBases, (byte) 126);
        this.cigarOperator = cigarOperator;
    }

    @Override
    public void setPixelRange(int s, int e) {
        this.pixelStart = s;
        this.pixelEnd = e;
    }

    @Override
    public boolean containsPixel(int x) {
        return x >= this.pixelStart && x <= this.pixelEnd;
    }

    @Override
    public int getPadding() {
        return padding;
    }

    @Override
    public char getCigarOperator() {
        return cigarOperator;
    }

    @Override
    public boolean contains(int position) {
        int offset = position - start;
        return offset >= 0 && offset < getLength();
    }

    @Override
    public int getBasesLength() {
        return basesLength;
    }

    @Override
    public int getBasesOffset() {
        return offset;
    }

    @Override
    public int getLength() {
        return basesLength + padding;
    }

    @Override
    public byte getBase(int offset) {
        return bases != null && offset < bases.length ? bases.getByte(offset) : 0;
    }

    /**
     * Return AlignmentBlock bases.  May be null
     *
     * @return
     */
    @Override
    public ByteSubarray getBases() {
        return bases == null ? EMPTY_ARRAY : bases;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public byte getQuality(int offset) {
        return (qualities == null || offset >= qualities.length) ? (byte) 126 : qualities.getByte(offset);
    }

    @Override
    public ByteSubarray getQualities() {
        return qualities == null ? EMPTY_ARRAY : qualities;
    }

    @Override
    public int getEnd() {
        return start + getLength();
    }

    @Override
    public boolean isSoftClip() {
        return softClipped;
    }

    public void setSoftClipped(boolean softClipped) {
        this.softClipped = softClipped;
    }


    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("[block ");
        sb.append(isSoftClip() ? "softClipped " : " ");
        sb.append(getStart());
        sb.append("-");
        sb.append(getEnd());
        sb.append(" ");
        if (bases != null) {
            sb.append(bases.getString());
        }
        sb.append("]");
        return sb.toString();
    }


    /**
     * Whether this AlignmentBlock has non-null bases
     *
     * @return
     */
    @Override
    public boolean hasBases() {
        return this.bases != null && this.bases.length > 0;
    }


    public void setPadding(int padding) {
        this.padding = padding;
    }
}
