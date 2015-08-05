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

import java.util.Arrays;

public class AlignmentBlock {

    private String chr;
    private int start;
    private byte[] bases;
    private int length = -1;
    public byte[] qualities;
    protected short[] counts;

    private boolean softClipped = false;

    private FlowSignalContext fContext = null;

    Alignment alignment;
    int offset;
    int end;

    /**
     * The reference genome we store mismatches to
     */
    private Genome genome;

    public AlignmentBlock(String chr, int start, byte[] bases, byte[] qualities) {
        this.chr = chr;
        this.start = start;
        this.bases = bases;
        this.length = bases.length;
        if (qualities == null || qualities.length < bases.length) {
            this.qualities = new byte[bases.length];
            Arrays.fill(this.qualities, (byte) 126);
        } else {
            this.qualities = qualities;
        }
        this.counts = null;
    }

    protected AlignmentBlock(String chr, int start, byte[] bases, byte[] qualities, FlowSignalContext fContext) {
        this(chr, start, bases, qualities);
        if (fContext != null && fContext.getNrSignals() == bases.length) {
            this.fContext = fContext;
        }
    }


    public boolean contains(int position) {
        int offset = position - start;
        return offset >= 0 && offset < getLength();
    }

    /**
     * Return AlignmentBlock bases.
     * May be null, which indicates they match the reference
     * @return
     */
    public byte[] getBases() {
        if(bases != null){
            return bases;
        }else{
            return getReferenceSequence();
        }
    }

    private byte[] getReferenceSequence() {
        return genome.getSequence(this.chr, getStart(), getEnd());
    }

    public int getLength() {
        return length;
    }

    public byte getBase(int offset) {
        return getBases()[offset];
    }

    public int getStart() {
        return start;
    }

    public byte getQuality(int offset) {
        return getQualities()[offset];

    }

    public byte[] getQualities() {
        return qualities;
    }

    public short getCount(int i) {
        return counts[i];
    }

    public void setCounts(short[] counts) {
        this.counts = counts;
    }

    public int getEnd() {
        return start + getLength();
    }

    public boolean isSoftClipped() {
        return softClipped;
    }

    public void setSoftClipped(boolean softClipped) {
        this.softClipped = softClipped;
    }

    public boolean hasCounts() {
        return counts != null;
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
     * Reduce so that we only store the mismatches between this block and reference
     * This may do nothing, if there are any mismatches we keep the original.
     * Note that we require an EXACT match, meaning ambiguity codes need to match
     * exactly. So if reference = 'N' and read = 'A', the full read sequence is stored.
     * @param genome
     */
    public void reduce(Genome genome){
        this.genome = genome;
        byte[] refBases = genome.getSequence(this.chr, getStart(), getEnd());
        //null refBases mostly happens in testing, but if we have no reference can't create mismatch
        if(refBases != null && this.bases != null){
            boolean match = false;
            for(int idx = 0; idx < refBases.length; idx++){
                match = refBases[idx] == this.bases[idx];
                if(!match) break;
            }
            if(match) this.bases = null;
        }
    }

    /**
     * Whether this AlignmentBlock has non-null bases
     * @return
     */
    public boolean hasBases() {
        return this.bases != null;
    }


    public FlowSignalSubContext getFlowSignalSubContext(int offset) {

        return  this.fContext == null ? null :
                new FlowSignalSubContext(this.fContext.getSignalForOffset(offset),
                        this.fContext.getBasesForOffset(offset), this.fContext.getFlowOrderIndexForOffset(offset));
    }


    public boolean hasFlowSignals() {
        return (null != this.fContext);
    }


}
