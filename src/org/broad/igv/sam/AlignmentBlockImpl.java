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

import org.broad.igv.feature.genome.Genome;

public class AlignmentBlockImpl implements AlignmentBlock {

    private String chr;
    private int start;
    private byte[] bases;
    private int length = -1;
    public byte[] qualities;

    private boolean softClipped = false;

    private FlowSignalContext fContext = null;

    Alignment alignment;
    int offset;
    int end;


    public AlignmentBlockImpl(String chr, int start, byte[] bases, byte[] qualities) {
        this.chr = chr;
        this.start = start;
        this.bases = bases;
        this.length = bases.length;
        this.qualities = qualities;
    }

    protected AlignmentBlockImpl(String chr, int start, byte[] bases, byte[] qualities, FlowSignalContext fContext) {
        this(chr, start, bases, qualities);
        if (fContext != null && fContext.getNrSignals() == bases.length) {
            this.fContext = fContext;
        }
    }

    @Override
    public boolean contains(int position) {
        int offset = position - start;
        return offset >= 0 && offset < getLength();
    }

    @Override
    public int getLength() {
        return length;
    }

    @Override
    public byte getBase(int offset) {
        return bases != null && offset < bases.length ? bases[offset] : 0;
    }

    /**
     * Return AlignmentBlock bases.  May be null
     *
     * @return
     */
    @Override
    public byte[] getBases() {
        return bases;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public byte getQuality(int offset) {
        return qualities == null || offset >= qualities.length ? (byte) 126 : qualities[offset];

    }

    @Override
    public byte[] getQualities() {
        return qualities;
    }

    @Override
    public int getEnd() {
        return start + getLength();
    }

    @Override
    public boolean isSoftClipped() {
        return softClipped;
    }

    @Override
    public void reduce(Genome genome) {

    }

    public void setSoftClipped(boolean softClipped) {
        this.softClipped = softClipped;
    }


    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("[block ");
        sb.append(isSoftClipped() ? "softClipped " : " ");
        sb.append(getStart());
        sb.append("-");
        sb.append(getEnd());
        sb.append(" ");
        for (int i = 0; i < bases.length; i++) {
            sb.append((char) bases[i]);
        }
        sb.append("]");
        return sb.toString();
    }


    /**
     * Whether this AlignmentBlock has non-null bases
     * @return
     */
    @Override
    public boolean hasBases() {
        return this.bases != null;
    }


    @Override
    public FlowSignalSubContext getFlowSignalSubContext(int offset) {

        return  this.fContext == null ? null :
                new FlowSignalSubContext(this.fContext.getSignalForOffset(offset),
                        this.fContext.getBasesForOffset(offset), this.fContext.getFlowOrderIndexForOffset(offset));
    }


    @Override
    public boolean hasFlowSignals() {
        return (null != this.fContext);
    }


}
