/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.apache.commons.lang.ArrayUtils;

import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AlignmentBlock {

    protected int start;
    protected byte[] bases;
    public byte[] qualities;
    protected short[] counts;

    private boolean softClipped = false;

    private SoftReference<byte[]> softBases;
    private SoftReference<byte[]> softQualities;
    private byte[] reference;

    /**
     * We save space by only storing the mismatches to the reference
     */
    private List<MismatchBlock> mismatches;

    public static AlignmentBlock getInstance(int start, byte[] bases, byte[] qualities) {

        return new AlignmentBlock(start, bases, qualities);
    }

    public static AlignmentBlock getInstance(int start, byte[] bases, byte[] qualities, FlowSignalContext fContext) {
        return new AlignmentBlockFS(start, bases, qualities, fContext);
    }

    protected AlignmentBlock(int start, byte[] bases, byte[] qualities) {
        this.start = start;
        this.bases = bases;
        if (qualities == null || qualities.length < bases.length) {
            this.qualities = new byte[bases.length];
            Arrays.fill(this.qualities, (byte) 126);
        } else {
            this.qualities = qualities;
        }
        this.counts = null;
    }

    public boolean contains(int position) {
        int offset = position - start;
        return offset >= 0 && offset < getLength();
    }

    public byte[] getBases() {
        if(bases != null) return bases;
        byte[] sbases = softBases.get();
        if(sbases == null){
            sbases = Arrays.copyOf(reference, reference.length);
            for(MismatchBlock mismatchBlock: this.mismatches){
                System.arraycopy(mismatchBlock.bases, 0, sbases, mismatchBlock.start - start, mismatchBlock.bases.length);
            }
            softBases = new SoftReference<byte[]>(sbases);
        }
        return sbases;
    }

    public int getLength() {
        if(bases != null) {
            return bases.length;
        }else if(reference != null) {
            return reference.length;
        }else{
            throw new IllegalStateException("Have no bases or reference");
        }
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
        if(qualities != null) return qualities;
        byte[] squals = softQualities.get();
        if(squals == null){
            squals = new byte[getLength()];
            Arrays.fill(squals, (byte) 126);
            for(MismatchBlock mismatchBlock: this.mismatches){
                System.arraycopy(mismatchBlock.qualities, 0, squals, mismatchBlock.start - start, mismatchBlock.qualities.length);
            }
            softQualities = new SoftReference<byte[]>(squals);
        }
        return squals;
    }

    public short[] getCounts() {
        return counts;
    }

    public short getCount(int i) {
        return counts[i];
    }

    public void setCounts(short[] counts) {
        this.counts = counts;
    }

    public int getEnd() {
        return start + getBases().length;
    }

    public boolean isSoftClipped() {
        return softClipped;
    }

    public void setSoftClipped(boolean softClipped) {
        this.softClipped = softClipped;
    }

    public boolean hasFlowSignals() {
        return false;
    }

    public boolean hasCounts() {
        return counts != null;
    }

    // Default implementation -- to be overridden
    public FlowSignalSubContext getFlowSignalSubContext(int offset) {
        return null;
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

    static List<MismatchBlock> createMismatchBlocks(int start, byte[] refBases, byte[] readBases, byte[] readQualities){
        List<MismatchBlock> mismatchBlocks = new ArrayList<MismatchBlock>();
        List<Byte> mismatches = null;
        List<Byte> qualities = null;
        assert readBases.length == refBases.length;
        int lastMMBlockStart = -1;
        for(int ii = 0; ii <= readBases.length; ii++){

            byte readBase = -1;
            byte refBase = -1;
            boolean atEnd = false;
            if(ii < readBases.length){
                readBase = readBases[ii];
                refBase = refBases[ii];
            }else{
                atEnd = true;
            }

            if(atEnd || AlignmentUtils.compareBases(refBase, readBase)){
                //Finish off last mismatch
                if(mismatches != null){
                    byte[] seq = ArrayUtils.toPrimitive(mismatches.toArray(new Byte[0]));
                    byte[] quals = ArrayUtils.toPrimitive(qualities.toArray(new Byte[0]));
                    MismatchBlock curMMBlock = new MismatchBlock(lastMMBlockStart, seq, quals);
                    mismatchBlocks.add(curMMBlock);
                    mismatches = null;
                    qualities = null;
                    lastMMBlockStart = -1;
                }
            }else{
                if(mismatches == null){
                    lastMMBlockStart = start + ii;
                    mismatches = new ArrayList<Byte>();
                    qualities = new ArrayList<Byte>();
                }
                mismatches.add(readBase);
                qualities.add(readQualities[ii]);
            }
        }

        return mismatchBlocks;
    }

    public List<MismatchBlock> getMismatches() {
        return mismatches;
    }

    /**
     * Reduce so that we only store the mismatches between this block and reference
     * @param refBases
     */
    public void reduce(byte[] refBases){
        mismatches = AlignmentBlock.createMismatchBlocks(getStart(), refBases, bases, qualities);
        this.reference = refBases;
        this.softBases = new SoftReference<byte[]>(this.bases);
        this.softQualities = new SoftReference<byte[]>(this.qualities);
        this.bases = null;
        qualities = null;
        counts = null;
    }


    public static class MismatchBlock{

        public final int start;
        public final byte[] bases;
        public final byte[] qualities;

        public MismatchBlock(int start, byte[] bases, byte[] qualities){
            this.start = start;
            this.bases = bases;
            this.qualities = qualities;
        }

    }
}
