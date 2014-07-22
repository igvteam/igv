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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author jrobinso
 */
public abstract class SAMAlignment implements Alignment {

    private static Logger log = Logger.getLogger(SAMAlignment.class);

    public static final char DELETE_CHAR = '-';
    public static final char SKIP_CHAR = '=';
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
    public static final String REDUCE_READS_TAG = "RR";

    private static final int READ_PAIRED_FLAG = 0x1;
    private static final int PROPER_PAIR_FLAG = 0x2;
    private static final int READ_UNMAPPED_FLAG = 0x4;
    protected static final int MATE_UNMAPPED_FLAG = 0x8;
    private static final int READ_STRAND_FLAG = 0x10;
    protected static final int MATE_STRAND_FLAG = 0x20;
    private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;
    private static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    private static final int DUPLICATE_READ_FLAG = 0x400;
    private static final int SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;
    /**
     * Converts a DNA integer value to its reverse compliment integer value.
     */
    protected static final char[] NT2COMP = {
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
            'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
    };
    private static final String FLOW_SIGNAL_TAG = "ZF";

    protected int flags;

    protected int start;  // <= Might differ from alignment start if soft clipping is considered
    protected int end;    // ditto
    protected int alignmentStart;
    protected int alignmentEnd;
    protected Color color = null;

    protected String readGroup;
    protected String library;
    protected String sample;

    String chr;
    int inferredInsertSize;
    int mappingQuality = 255;  // 255 by default
    ReadMate mate;
    String readName;
    AlignmentBlock[] alignmentBlocks;
    AlignmentBlock[] insertions;
    char[] gapTypes;
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();

    protected String cigarString;
    protected String readSequence;

    protected String mateSequence = null;
    protected String pairOrientation = "";
    private Strand firstOfPairStrand;
    private Strand secondOfPairStrand;

    public SAMAlignment() {
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

    public Color getColor() {
        return color;
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

    public boolean isNegativeStrand() {
        return (flags & READ_STRAND_FLAG) != 0;
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
    /**
     * Set pair strands.  Used for strand specific libraries to recover strand of
     * originating fragment.
     */
    protected void setPairStrands() {

        if (isPaired()) {
            if (isFirstOfPair()) {
                firstOfPairStrand = getReadStrand();
            } else {
                // If we have a mate, the mate must be the firstOfPair
                ReadMate mate = getMate();
                if (mate != null && mate.isMapped()) {
                    firstOfPairStrand = mate.getStrand();
                } else {
                    // No Mate, or mate is not mapped, FOP strand is not defined
                    firstOfPairStrand = Strand.NONE;
                }
            }

            if (isSecondOfPair()) {
                secondOfPairStrand = isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
            } else {
                ReadMate mate = getMate();
                if (mate.isMapped() && isProperPair()) {
                    secondOfPairStrand = mate.isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
                } else {
                    // No Mate, or mate is not mapped, FOP strand is not defined
                    secondOfPairStrand = Strand.NONE;
                }
            }

        } else {
            // This alignment is not paired -- by definition "firstOfPair" is this alignment
            firstOfPairStrand = getReadStrand();
            secondOfPairStrand = Strand.NONE;
        }
    }


    private static boolean operatorIsMatch(boolean showSoftClipped, char operator) {
        return operator == MATCH || operator == PERFECT_MATCH || operator == MISMATCH
                || (showSoftClipped && operator == SOFT_CLIP);
    }



    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     *
     * @param cigarString
     * @param readBases
     * @param readBaseQualities
     * @param flowSignals              from the FZ tag, null if not present
     * @param flowOrder                from the RG.FO header tag, null if not present
     * @param flowOrderStart
     */
    protected void createAlignmentBlocks(String cigarString, byte[] readBases, byte[] readBaseQualities,
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
                            getChr(), blockStart, fromIdx, op.nBases, true);

                    if (op.operator == SOFT_CLIP) {
                        block.setSoftClipped(true);
                    }
                    alignmentBlocks[blockIdx++] = block;

                    fromIdx += op.nBases;
                    blockStart += op.nBases;

                    if (operatorIsMatch(showSoftClipped, prevOp)) {
                        gapTypes[gapIdx++] = ZERO_GAP;
                    }

                } else if (op.operator == DELETION || op.operator == SKIPPED_REGION) {
                    blockStart += op.nBases;
                    gapTypes[gapIdx++] = op.operator;
                } else if (op.operator == INSERTION) {
                    // This gap is between blocks split by insertion.   It is a zero
                    // length gap but must be accounted for.
                    gapTypes[gapIdx++] = ZERO_GAP;
                    AlignmentBlock block = buildAlignmentBlock(fBlockBuilder, readBases, readBaseQualities,
                            getChr(), blockStart, fromIdx, op.nBases, false);

                    insertions[insertionIdx++] = block;
                    fromIdx += op.nBases;
                } else if (op.operator == PADDING) {
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
                                                      byte[] readBaseQualities, String chr, int blockStart,
                                                      int fromIdx, int nBases, boolean checkNBasesAvailable) {

        byte[] blockBases = new byte[nBases];
        byte[] blockQualities = new byte[nBases];
        short[] blockCounts = new short[nBases];

        // TODO -- represent missing sequence ("*") explicitly for efficiency.
        int nBasesAvailable = nBases;
        if (checkNBasesAvailable) {
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
        if (checkNBasesAvailable) {
            nBasesAvailable = readBaseQualities.length - fromIdx;
        }
        if (readBaseQualities == null || readBaseQualities.length == 0 || nBasesAvailable < nBases) {
            Arrays.fill(blockQualities, (byte) 126);
        } else {
            System.arraycopy(readBaseQualities, fromIdx, blockQualities, 0, nBases);
        }

        AlignmentBlock block;
        if (fBlockBuilder != null) {
            block = AlignmentBlock.getInstance(chr, blockStart, blockBases, blockQualities,
                    fBlockBuilder.getFlowSignalContext(readBases, fromIdx, nBases));
        } else {
            block = AlignmentBlock.getInstance(chr, blockStart, blockBases, blockQualities);
        }

        return block;
    }


    private static void bufAppendFlowSignals(AlignmentBlock block, StringBuffer buf, int offset) {
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

    public String getClipboardString(double location) {
        return getValueString(location, null);
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

    public boolean isFirstOfPair() {
        return  isPaired() && (flags & FIRST_OF_PAIR_FLAG) != 0;
    }

    public boolean isSecondOfPair() {
        return isPaired() && (flags & SECOND_OF_PAIR_FLAG) != 0;
    }

    /**
     * @return the unclippedStart
     */
    public int getAlignmentStart() {
        return alignmentStart;
    }

    public String getCigarString() {
        return cigarString;
    }
    public int getReadLength() {
        return readSequence.length();
    }

    public boolean isDuplicate() {
        return (flags & DUPLICATE_READ_FLAG) != 0;
    }

    public boolean isMapped() {
        return (flags & READ_UNMAPPED_FLAG) == 0;
    }


    public boolean isPaired() {
        return (flags & READ_PAIRED_FLAG) != 0;
    }

    public boolean isProperPair(){
      return  ((flags & READ_PAIRED_FLAG) != 0) && ((flags & PROPER_PAIR_FLAG) != 0);
    }

    public boolean isSmallInsert() {
        int absISize = Math.abs(getInferredInsertSize());
        return absISize > 0 && absISize <= getReadLength();
    }

    public float getScore() {
        return getMappingQuality();
    }

    /**
     * @param mate the mate to set
     */
    public void setMate(ReadMate mate) {
        this.mate = mate;
    }

    @Override
    public boolean isSupplementary() {
        return (flags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }

    /**
     * @return the alignmentEnd
     */
    public int getAlignmentEnd() {
        return alignmentEnd;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getSample() {
        return sample;
    }

    public String getReadGroup() {
        return readGroup;
    }


    public abstract Object getAttribute(String key);

    public String getLibrary() {
        return library;
    }

    @Override
    public char[] getGapTypes() {
        return gapTypes;
    }


    @Override
    public void finish() {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        for (AlignmentBlock block : alignmentBlocks) {
            block.reduce(genome);
        }
    }


    public boolean isVendorFailedRead() {
        return (flags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public String getPairOrientation() {
        return pairOrientation;
    }


    public String getReadSequence() {
        return readSequence;
    }

    /**
     * Use blocks to recreate read sequence.
     * As of this comment writing, we don't keep a block
     * for hard-clipped bases, so this won't match what's in the file
     *
     * @return
     */
    String buildReadSequenceFromBlocks() {
        String readSeq = "";
        for (AlignmentBlock block : getAlignmentBlocks()) {
            readSeq += new String(block.getBases());
        }
        return readSeq;
    }


    @Override
    public boolean isPrimary() {
        return (flags & NOT_PRIMARY_ALIGNMENT_FLAG) == 0;
    }



    @Override
    public void setMateSequence(String sequence) {
        this.mateSequence = sequence;
    }

    /**
     * Return the strand of the read marked "first-in-pair" for a paired alignment. This method can return
     * Strand.NONE if the end marked first is unmapped.
     *
     * @return strand of first-of-pair
     */
    public Strand getFirstOfPairStrand() {
        return firstOfPairStrand;
    }

    /**
     * Return the strand of the read marked "second-in-pair" for a paired alignment.  The strand is
     * undefined (Strand.NONE) for non-paired alignments
     *
     * @return strand of second-of-pair
     */
    public Strand getSecondOfPairStrand() {
        return secondOfPairStrand;
    }

    /**
     * @return start index in the flow signal as specified by the ZF tag, or -1 if not present
     * or non-numeric
     */
    public int getFlowSignalsStart() {
        Object attribute = getAttribute(FLOW_SIGNAL_TAG); // NB: from a TMAP optional tag
        int toRet = -1;
        if (attribute != null && attribute instanceof Integer) {
            toRet = (Integer) attribute;
        }
        return toRet;
    }


    protected void setPairOrientation() {

        /*
                if (record.getReadPairedFlag() &&
                !record.getReadUnmappedFlag() &&
                !record.getMateUnmappedFlag() &&
                record.getReferenceName().equals(record.getMateReferenceName())) {

         */

        if (isPaired() && isMapped() && mate.isMapped() && getChr().equals(mate.getChr())) {   // && name === mate.name

            char s1 = isNegativeStrand() ? 'R' : 'F';
            char s2 = mate.isNegativeStrand() ? 'R' : 'F';
            char o1 = ' ';
            char o2 = ' ';
            if (isFirstOfPair() ) {
                o1 = '1';
                o2 = '2';
            } else if (isSecondOfPair()) {
                o1 = '2';
                o2 = '1';
            }

            final char[] tmp = new char[4];
            int isize = inferredInsertSize;
            int estReadLen = getAlignmentEnd() - getAlignmentStart();
            if (isize == 0) {
                //isize not recorded.  Need to estimate.  This calculation was validated against an Illumina
                // -> <- library bam.
                int estMateEnd = getAlignmentStart() < mate.getStart() ?
                        getMate().getStart() + estReadLen : mate.getStart() - estReadLen;
                isize = estMateEnd - getAlignmentStart();
            }

            //if (isize > estReadLen) {
            if (isize > 0) {
                tmp[0] = s1;
                tmp[1] = o1;
                tmp[2] = s2;
                tmp[3] = o2;

            } else {
                tmp[2] = s1;
                tmp[3] = o1;
                tmp[0] = s2;
                tmp[1] = o2;
            }
            // }
            pairOrientation = new String(tmp);
        }
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
