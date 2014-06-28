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

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author jrobinso
 */
public abstract class AbstractAlignment implements Alignment {

    public static final char DELETE_CHAR = '-';
    public static final char SKIP_CHAR = '=';
    public static final String REDUCE_READS_TAG = "RR";
    protected static final String FLOW_SIGNAL_TAG = "ZF";
    private static Logger log = Logger.getLogger(AbstractAlignment.class);

    public static final char MATCH = 'M';
    public static final char PERFECT_MATCH = '=';
    public static final char MISMATCH = 'X';
    public static final char INSERTION = 'I';
    public static final char DELETION = 'D';
    public static final char SKIPPED_REGION = 'N';
    public static final char SOFT_CLIP = 'S';
    public static final char HARD_CLIP = 'H';
    public static final char PADDING = 'P';
    public static final char ZERO_GAP = 'O';
    protected int start;  // <= Might differ from alignment start if soft clipping is considered
    protected int end;    // ditto
    String chr;
    int inferredInsertSize;
    int mappingQuality = 255;  // 255 by default
    ReadMate mate;
    String readName;
    AlignmentBlock[] alignmentBlocks;
    AlignmentBlock[] insertions;
    char[] gapTypes;
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();
    public boolean negativeStrand;
    public boolean duplicate;
    public boolean mapped;
    public int readLength;
    public boolean paired;
    public boolean properPair;
    public boolean firstOfPair;
    public boolean secondOfPair;
    public String cigarString;
    public String readSequence;
    public boolean primary;
    public boolean supplementary;

    public AbstractAlignment() {
    }

    public String getChromosome() {
        return getChr();
    }

    public String getChr() {
        return chr;
    }

    public String getDescription() {
        return getReadName();
    }

    public ReadMate getMate() {
        return mate;
    }

    public String getMateSequence() {
        return null;
    }

    public String getReadName() {
        return readName;
    }

    public int getMappingQuality() {
        return mappingQuality;
    }

    public int getInferredInsertSize() {
        return inferredInsertSize;
    }

    public AlignmentBlock[] getAlignmentBlocks() {
        return alignmentBlocks;
    }

    public AlignmentBlock[] getInsertions() {
        return insertions;
    }

    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     *
     * @param cigarString
     * @param readBases
     * @param readBaseQualities
     * @param readRepresentativeCounts the representative counts of each base in the read (translated from the reduce reads tag)
     * @param flowSignals              from the FZ tag, null if not present
     * @param flowOrder                from the RG.FO header tag, null if not present
     * @param flowOrderStart
     */
    protected void createAlignmentBlocks(String cigarString, byte[] readBases, byte[] readBaseQualities, short[] readRepresentativeCounts,
                                         short[] flowSignals, String flowOrder, int flowOrderStart) {

        boolean showSoftClipped = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED);

        int nInsertions = 0;
        int nBlocks = 0;

        java.util.List<CigarOperator> operators = new ArrayList();
        StringBuffer buffer = new StringBuffer(4);

        if (cigarString.equals("*")) {
            alignmentBlocks = new AlignmentBlock[1];
            alignmentBlocks[0] = new AlignmentBlock(getChr(), getStart(), readBases, readBaseQualities);
            return;
        }

        // Create list of cigar operators
        boolean firstOperator = true;
        int softClippedBaseCount = 0;
        int nGaps = 0;
        char prevOp = 0;
        for (int i = 0; i < cigarString.length(); i++) {
            char next = cigarString.charAt(i);
            if (Character.isDigit(next)) {
                buffer.append(next);
            } else {
                char op = next;
                if (op == HARD_CLIP) {
                    buffer = new StringBuffer(4);
                    continue;  // Just skip hardclips
                }
                int nBases = Integer.parseInt(buffer.toString());
                if (operatorIsMatch(showSoftClipped, op)) {
                    if (operatorIsMatch(showSoftClipped, prevOp)) {
                        nGaps++;   // Consecutive Ms
                    }
                    nBlocks++;

                } else if (op == DELETION || op == SKIPPED_REGION) {
                    nGaps++;
                } else if (op == INSERTION) {
                    nInsertions++;
                    nGaps++; // "virtual" gap, account for artificial block split @ insertion
                } else if (op == PADDING) {
                    nGaps++;
                }


                if (firstOperator && op == SOFT_CLIP) {
                    softClippedBaseCount += nBases;
                }

                operators.add(new CigarOperator(nBases, op));
                buffer = new StringBuffer(4);

                prevOp = op;
                firstOperator = false;
            }
        }

        alignmentBlocks = new AlignmentBlock[nBlocks];
        insertions = new AlignmentBlock[nInsertions];
        if (nGaps > 0) {
            gapTypes = new char[nGaps];
        }

        // Adjust start to include soft clipped bases a
        if (showSoftClipped) {
            start -= softClippedBaseCount;
        }
        int fromIdx = showSoftClipped ? 0 : softClippedBaseCount;
        int blockStart = start;

        // Create blocks
        int blockIdx = 0;
        int insertionIdx = 0;
        int gapIdx = 0;
        FlowSignalContextBuilder fBlockBuilder = null;
        if (null != flowSignals) {
            if (0 < readBases.length) {
                fBlockBuilder = new FlowSignalContextBuilder(flowSignals, flowOrder, flowOrderStart, readBases, fromIdx, this.isNegativeStrand());
            }
        }
        prevOp = 0;
        for (CigarOperator op : operators) {
            try {

                if (op.operator == HARD_CLIP) {
                    continue;
                }
                if (operatorIsMatch(showSoftClipped, op.operator)) {

                    AlignmentBlock block = buildAlignmentBlock(fBlockBuilder, readBases, readBaseQualities,
                            readRepresentativeCounts, getChr(), blockStart, fromIdx, op.nBases, true);

                    if (op.operator == SOFT_CLIP) {
                        block.setSoftClipped(true);
                    }
                    alignmentBlocks[blockIdx++] = block;

                    fromIdx += op.nBases;
                    blockStart += op.nBases;

                    if (operatorIsMatch(showSoftClipped, prevOp)) {
                        gapTypes[gapIdx++] = ZERO_GAP;
                    }

                } else if (op.operator == DELETION || op.operator == SKIPPED_REGION ) {
                    blockStart += op.nBases;
                    gapTypes[gapIdx++] = op.operator;
                } else if (op.operator == INSERTION) {
                    // This gap is between blocks split by insertion.   It is a zero
                    // length gap but must be accounted for.
                    gapTypes[gapIdx++] = ZERO_GAP;
                    AlignmentBlock block = buildAlignmentBlock(fBlockBuilder, readBases, readBaseQualities,
                            readRepresentativeCounts, getChr(), blockStart, fromIdx, op.nBases, false);

                    insertions[insertionIdx++] = block;
                    fromIdx += op.nBases;
                } else if (op.operator == PADDING){
                    //Padding represents a deletion against the padded reference
                    //But we don't have the padded reference
                    gapTypes[gapIdx++] = ZERO_GAP;
                }
            } catch (Exception e) {
                log.error("Error processing CIGAR string", e);
            }
            prevOp = op.operator;
        }

        // Check for soft clipping at end
        if (showSoftClipped && operators.size() > 0) {
            CigarOperator last = operators.get(operators.size() - 1);
            if (last.operator == SOFT_CLIP) {
                end += last.nBases;
            }
        }
    }

    private static AlignmentBlock buildAlignmentBlock(FlowSignalContextBuilder fBlockBuilder, byte[] readBases,
                                               byte[] readBaseQualities,
                                               short[] readRepresentativeCounts, String chr, int blockStart,
                                               int fromIdx, int nBases, boolean checkNBasesAvailable){

        byte[] blockBases = new byte[nBases];
        byte[] blockQualities = new byte[nBases];
        short[] blockCounts = new short[nBases];

        // TODO -- represent missing sequence ("*") explicitly for efficiency.
        int nBasesAvailable = nBases;
        if(checkNBasesAvailable){
            nBasesAvailable = readBases.length - fromIdx;
        }
        if (readBases == null || readBases.length == 0) {
            Arrays.fill(blockBases, (byte) '=');
        } else if (nBasesAvailable < nBases) {
            Arrays.fill(blockBases, (byte) '?');
        } else {
            System.arraycopy(readBases, fromIdx, blockBases, 0, nBases);
        }

        nBasesAvailable = nBases;
        if(checkNBasesAvailable){
            nBasesAvailable = readBaseQualities.length - fromIdx;
        }
        if (readBaseQualities == null || readBaseQualities.length == 0 || nBasesAvailable < nBases) {
            Arrays.fill(blockQualities, (byte) 126);
        } else {
            System.arraycopy(readBaseQualities, fromIdx, blockQualities, 0, nBases);
        }

        if (readRepresentativeCounts != null) {
            System.arraycopy(readRepresentativeCounts, fromIdx, blockCounts, 0, nBases);
        }

        AlignmentBlock block;
        if (fBlockBuilder != null) {
            block = AlignmentBlock.getInstance(chr, blockStart, blockBases, blockQualities,
                    fBlockBuilder.getFlowSignalContext(readBases, fromIdx, nBases));
        } else {
            block = AlignmentBlock.getInstance(chr, blockStart, blockBases, blockQualities);
        }
        if (readRepresentativeCounts != null) {
            block.setCounts(blockCounts);
        }
        return block;
    }

    private static boolean operatorIsMatch(boolean showSoftClipped, char operator) {
        return operator == MATCH || operator == PERFECT_MATCH || operator == MISMATCH
                || (showSoftClipped && operator == SOFT_CLIP);
    }

    public boolean isNegativeStrand() {
        return negativeStrand;
    }

    public void setNegativeStrand(boolean negativeStrand) {
        this.negativeStrand = negativeStrand;
    }

    public boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    public byte getBase(double position) {
        int basePosition = (int) position;
        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte base = block.getBases()[offset];
                return base;
            }
        }
        return 0;
    }

    public byte getPhred(double position) {
        int basePosition = (int) position;
        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte qual = block.getQuality(offset);
                return qual;
            }
        }
        return 0;
    }

    private void bufAppendFlowSignals(AlignmentBlock block, StringBuffer buf, int offset) {
        if (block.hasFlowSignals()) {
            // flow signals
            int i, j, n = 0;
            FlowSignalSubContext f = block.getFlowSignalSubContext(offset);
            if (null != f && null != f.getSignals() && null != f.getBases()) {
                buf.append("FZ = ");
                StringBuffer spos = new StringBuffer();
                spos.append("Flow position = ").append(f.getFlowOrderIndex());

                for (i = 0; i < f.getNrSignalTypes(); i++) {
                    short[] signals = f.getSignalsOfType(i);
                    char[] bases = f.getBasesOfType(i);
                    if (null != signals && 0 < signals.length) {
                        if (1 == i) {
                            if (0 < n) {
                                buf.append(",");
                            }
                            buf.append("[");
                        }
                        for (j = 0; j < signals.length; j++) {
                            if (1 != i && 0 < n) {
                                buf.append(",");
                            }
                            buf.append(bases[j]);
                            buf.append(signals[j]);

                            n++;
                        }
                        if (1 == i) {
                            buf.append("]");
                        }
                    }
                }
                buf.append("<br>").append(spos);
                buf.append("<br>");
                // maybe also add flow order?                
            }
        }
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = null;

        // First check insertions.  Position is zero based, block coords 1 based
        if (this.insertions != null) {
            for (AlignmentBlock block : this.insertions) {
                double insertionLeft = block.getStart() - .25;
                double insertionRight = block.getStart() + .25;
                if (position > insertionLeft && position < insertionRight) {
                    if (block.hasFlowSignals()) {
                        int offset;
                        buf = new StringBuffer();
                        buf.append("Insertion: " + new String(block.getBases()) + "<br>");
                        buf.append("Base phred quality = ");
                        for (offset = 0; offset < block.getLength(); offset++) {
                            byte quality = block.getQuality(offset);
                            if (0 < offset) {
                                buf.append(",");
                            }
                            buf.append(quality);
                        }
                        buf.append("<br>");
                        for (offset = 0; offset < block.getLength(); offset++) {
                            byte base = block.getBase(offset);
                            buf.append((char) base + ": ");
                            bufAppendFlowSignals(block, buf, offset);
                        }
                        buf.append("----------------------"); // NB: no <br> required
                        return buf.toString();
                    } else {
                        return "Insertion: " + new String(block.getBases());
                    }
                }
            }
        }

        buf = new StringBuffer();

        String sample = getSample();
        if (sample != null) {
            buf.append("Sample = " + sample + "<br>");
        }
        String readGroup = getReadGroup();
        if (sample != null) {
            buf.append("Read group = " + readGroup + "<br>");
        }
        buf.append("----------------------" + "<br>");

        int basePosition = (int) position;
        buf.append("Read name = " + getReadName() + "<br>");
        buf.append("Location = " + getChr() + ":" + DECIMAL_FORMAT.format(1 + (long) position) + "<br>");
        buf.append("Alignment start = " + DECIMAL_FORMAT.format(getAlignmentStart() + 1) + " (" + (isNegativeStrand() ? "-" : "+") + ")<br>");
        buf.append("Cigar = " + getCigarString() + "<br>");
        buf.append("Mapped = " + (isMapped() ? "yes" : "no") + "<br>");
        buf.append("Mapping quality = " + getMappingQuality() + "<br>");
        buf.append("----------------------" + "<br>");

        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte base = block.getBase(offset);
                byte quality = block.getQuality(offset);
                buf.append("Base = " + (char) base + "<br>");
                buf.append("Base phred quality = " + quality + "<br>");
                if (block.hasCounts()) {
                    buf.append("Count = " + block.getCount(offset) + "<br>");
                }
                // flow signals
                if (block.hasFlowSignals()) {
                    bufAppendFlowSignals(block, buf, offset);
                }
            }
        }

        if (this.isPaired()) {
            buf.append("----------------------" + "<br>");
            buf.append("Mate start = " + getMate().positionString() + "<br>");
            buf.append("Mate is mapped = " + (getMate().isMapped() ? "yes" : "no") + "<br>");
            //buf.append("Pair is proper = " + (getProperPairFlag() ? "yes" : "no") + "<br>");
            if (getChr().equals(getMate().getChr())) {
                buf.append("Insert size = " + getInferredInsertSize() + "<br>");
            }
            if (getPairOrientation().length() > 0) {
                buf.append("Pair orientation = " + getPairOrientation() + "<br>");
            }
        }
        buf.append("----------------------");
        return buf.toString();
    }

    public abstract String getCigarString();

    public abstract boolean isMapped();

    public abstract boolean isPaired();

    public abstract boolean isProperPair();

    public boolean isSmallInsert() {
        int absISize = Math.abs(getInferredInsertSize());
        return absISize > 0 && absISize <= getReadLength();
    }

    public int getReadLength() {
        return getReadSequence().length();
    }

    public float getScore() {
        return getMappingQuality();
    }

    /**
     * @param mappingQuality the mappingQuality to set
     */
    public void setMappingQuality(int mappingQuality) {
        this.mappingQuality = mappingQuality;
    }

    /**
     * @param inferredInsertSize the inferredInsertSize to set
     */
    public void setInferredInsertSize(int inferredInsertSize) {
        this.inferredInsertSize = inferredInsertSize;
    }

    /**
     * @param mate the mate to set
     */
    public void setMate(ReadMate mate) {
        this.mate = mate;
    }

    public String getReadGroup() {
        return null;
    }

    public String getLibrary() {
        return null;
    }

    public String getClipboardString(double location) {
        return getValueString(location, null);
    }

    public char[] getGapTypes() {
        return null;
    }

    public Object getAttribute(String key) {
        return null;
    }

    public void setMateSequence(String sequence) {
        // ignore by default
    }

    public String getPairOrientation() {
        return "";
    }

    public boolean isVendorFailedRead() {
        return false;
    }

    public Color getColor() {
        return null;
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public void finish() {
        //default operation is nothing
    }

    @Override
    public boolean isPrimary() {
        return true;
    }

    @Override
    public boolean isSupplementary() {
        return false;
    }

    static class CigarOperator {

        int nBases;
        char operator;

        /**
         * Constructs ...
         *
         * @param nBases
         * @param operator
         */
        public CigarOperator(int nBases, char operator) {
            this.nBases = nBases;
            this.operator = operator;
        }
    }
}
