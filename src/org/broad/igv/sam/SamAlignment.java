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
package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author jrobinso
 */
public class SamAlignment extends AbstractAlignment implements Alignment {

    private static Logger log = Logger.getLogger(SamAlignment.class);

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
    private int start;  // <= Might differ from alignment start if soft clipping is considered
    private int end;    // ditto
    private int alignmentStart;
    private int alignmentEnd;
    private int readLength;
    private boolean readNegativeStrandFlag;
    private boolean duplicateReadFlag;
    private boolean readUnmappedFlag;
    private boolean readPairedFlag;
    private boolean properPairFlag;
    private SAMRecord record;
    private String cigarString;
    private String readSequence;
    private boolean firstRead = false;
    private boolean secondRead = false;
    private String mateSequence = null;
    private String pairOrientation = "";
    private Color defaultColor = AlignmentRenderer.grey1;
    private String readGroup;
    private String library;
    private String sample;

    private boolean firstInPair;
    private Strand firstOfPairStrand;
    private Strand secondOfPairStrand;

    /**
     * Converts a DNA integer value to its reverse compliment integer value.
     */
    protected static final char NT2COMP[] = {
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

    // Default constructor to support unit tests
    //SamAlignment() {}


    public SamAlignment(SAMRecord record) {
        String flowOrder = null;
        String keySequence = null;

        this.record = record;

        String refName = record.getReferenceName();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.chr = genome == null ? refName : genome.getChromosomeAlias(refName);

        // SAMRecord is 1 based inclusive.  IGV is 0 based exclusive.
        this.alignmentStart = record.getAlignmentStart() - 1;
        this.start = this.alignmentStart;   // might be modified later for soft clipping
        this.alignmentEnd = Math.max(alignmentStart, record.getAlignmentEnd());
        this.end = alignmentEnd;   // might be modified later for soft clipping
        this.cigarString = record.getCigarString();
        this.setMappingQuality(record.getMappingQuality());
        this.readName = record.getReadName().trim();
        this.readNegativeStrandFlag = record.getReadNegativeStrandFlag();
        this.duplicateReadFlag = record.getDuplicateReadFlag();
        this.readUnmappedFlag = record.getReadUnmappedFlag();
        this.readPairedFlag = record.getReadPairedFlag();
        this.setInferredInsertSize(record.getInferredInsertSize());
        this.readSequence = record.getReadString();
        this.readLength = record.getReadLength();
        this.firstInPair = record.getReadPairedFlag() ? record.getFirstOfPairFlag() : true;

        setMatePair(record, genome);
        setPairOrientation(record);
        setPairStrands();

        SAMFileHeader header = record.getHeader();
        if (header != null) {
            readGroup = (String) record.getAttribute("RG");
            if (readGroup != null) {
                SAMReadGroupRecord rgRec = header.getReadGroup(readGroup);
                if (rgRec != null) {
                    String platform = rgRec.getPlatform();
                    this.sample = rgRec.getSample();
                    this.library = rgRec.getLibrary();
                    flowOrder = rgRec.getFlowOrder();
                    keySequence = rgRec.getKeySequence();
                }
            }
        }

        createAlignmentBlocks(record.getCigarString(), record.getReadBases(), record.getBaseQualities(),
                getFlowSignals(record, flowOrder, keySequence), flowOrder, this.getFlowSignalsStart(record));

        Object colorTag = record.getAttribute("YC");
        if (colorTag != null) {
            try {
                defaultColor = ColorUtilities.stringToColor(colorTag.toString());
            } catch (Exception e) {
                log.error("Error interpreting color tag: " + colorTag, e);
                defaultColor = AlignmentRenderer.grey1;
            }
        }
    }      // End constructor


    private void setMatePair(SAMRecord record, Genome genome) {
        if (record.getReadPairedFlag()) {
            String mateReferenceName = record.getMateReferenceName();
            String mateChr = genome == null ? mateReferenceName : genome.getChromosomeAlias(mateReferenceName);
            this.properPairFlag = record.getProperPairFlag();
            this.setMate(new ReadMate(mateChr,
                    record.getMateAlignmentStart(),
                    record.getMateNegativeStrandFlag(),
                    record.getMateUnmappedFlag()));
            this.firstRead = record.getFirstOfPairFlag();
            this.secondRead = record.getSecondOfPairFlag();
        }

    }

    private void setPairOrientation(SAMRecord record) {
        if (record.getReadPairedFlag() &&
                !readUnmappedFlag &&
                !record.getMateUnmappedFlag() &&
                record.getReferenceName().equals(record.getMateReferenceName())) {

            char s1 = record.getReadNegativeStrandFlag() ? 'R' : 'F';
            char s2 = record.getMateNegativeStrandFlag() ? 'R' : 'F';
            char o1 = ' ';
            char o2 = ' ';
            if (record.getFirstOfPairFlag()) {
                o1 = '1';
                o2 = '2';
            } else if (record.getSecondOfPairFlag()) {
                o1 = '2';
                o2 = '1';
            }

            final char[] tmp = new char[4];
            if (record.getInferredInsertSize() > 0) {
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
            pairOrientation = new String(tmp);
        }
    }

    /**
     * Set pair strands.  Used for strand specific libraries to recover strand of
     * originating fragment.
     */
    private void setPairStrands() {

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

    // Helper class to get the flow signals context
    private class FlowOrderAndSignalsBlockHelper {

        private short[] flowSignals = null;
        private String flowOrder = null;
        private int flowSignalsIndex = -1;
        private int flowOrderIndex = -1;
        private int prevFlowSignalsStart = -1;
        private int prevFlowSignalsEnd = -1;
        private boolean readNegativeStrandFlag;
        private boolean[] incorporations = null; // required for the reverse strand

        public FlowOrderAndSignalsBlockHelper(short[] flowSignals, String flowOrder, int flowOrderStart, byte[] readBases, boolean readNegativeStrandFlag) {
            if (null == flowSignals || null == flowOrder || flowOrderStart < 0) {
                return;
            }

            this.flowSignals = flowSignals;
            this.flowOrder = flowOrder;
            this.flowOrderIndex = flowOrderStart;
            this.flowSignalsIndex = 0; // NB: the key sequence/barcode sequence should have been remove for the signals already
            this.readNegativeStrandFlag = readNegativeStrandFlag;

            // init
            if (this.readNegativeStrandFlag) {
                int i;
                this.incorporations = new boolean[this.flowSignals.length];
                // go to the end of the signals
                for (i = readBases.length - 1; 0 <= i; i--) {
                    while (this.flowOrder.charAt(this.flowOrderIndex) != NT2COMP[readBases[i]]) {
                        this.flowOrderIndex++;
                        this.flowSignalsIndex++;
                        this.incorporations[this.flowSignalsIndex] = false;
                        if (this.flowOrder.length() <= this.flowOrderIndex) {
                            this.flowOrderIndex = 0;
                        }
                    }
                    this.incorporations[this.flowSignalsIndex] = true;
                }
                this.prevFlowSignalsStart = this.flowSignalsIndex + 1;
                this.prevFlowSignalsEnd = this.flowSignals.length - 1;
            } else {
                this.prevFlowSignalsStart = this.prevFlowSignalsEnd = 0;
                while (this.flowOrder.charAt(this.flowOrderIndex) != readBases[0]) {
                    this.flowOrderIndex++;
                    this.flowSignalsIndex++;
                    if (this.flowOrder.length() <= this.flowOrderIndex) {
                        this.flowOrderIndex = 0;
                    }
                }
                this.prevFlowSignalsEnd = this.flowSignalsIndex - 1;
            }
        }

        // TODO:
        // - support IUPAC bases
        // - support lower/upper cases (is this necessary)?
        public short[][][] createBlockFlowSignals(byte[] readBases, int fromIdx, int nBases) {
            int i, idx;
            short[][][] blockFlowSignals = null;

            if (null == this.flowSignals) {
                return null;
            }

            blockFlowSignals = new short[nBases][][];
            //Default value
            Arrays.fill(blockFlowSignals, null);

            // NB: should be at the first base of a HP
            // Go through the bases
            i = fromIdx;
            idx = 0;
            while (0 <= this.flowSignalsIndex && this.flowSignalsIndex < this.flowSignals.length && i < fromIdx + nBases) {
                short s = this.flowSignals[this.flowSignalsIndex];
                int nextFlowSignalsStart = -1, nextFlowSignalsEnd = -1;
                int j = i + 1;
                if (j < readBases.length) {
                    if (this.readNegativeStrandFlag) {
                        nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                        // NB: loop condition is not symmetric to the forward, as we must respect the directionality of sequencing.
                        // For example, if our flow order is TACAG, and our read bases are TAG, then the flow signal vector is 
                        // approximately 100,100,0,0,100.  Since we move in the reverse direction with respect to the flow signal 
                        // vector we must pre-compute where the flows incorporations are expected to occur, instead of just looking 
                        // for the next flow that matches our next read base (we would place the A incorporation flow in the fourth flow,
                        // which is wrong).
                        while (!this.incorporations[this.flowSignalsIndex] ||
                                this.flowOrder.charAt(this.flowOrderIndex) != NT2COMP[readBases[j]]) { // NB: malicious input can cause infinite loops here
                            this.flowOrderIndex--;
                            this.flowSignalsIndex--;
                            if (this.flowOrderIndex < 0) {
                                this.flowOrderIndex = this.flowOrder.length() - 1;
                            }
                        }
                        nextFlowSignalsStart = this.flowSignalsIndex + 1;
                    } else {
                        nextFlowSignalsStart = this.flowSignalsIndex + 1;
                        while (this.flowOrder.charAt(this.flowOrderIndex) != readBases[j]) { // NB: malicious input can cause infinite loops here
                            this.flowOrderIndex++;
                            this.flowSignalsIndex++;
                            if (this.flowOrder.length() <= this.flowOrderIndex) {
                                this.flowOrderIndex = 0;
                            }
                        }
                        nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                    }
                }
                // set-up block
                blockFlowSignals[idx] = new short[3][];
                // this.previous context
                if (0 <= this.prevFlowSignalsStart && this.prevFlowSignalsStart <= this.prevFlowSignalsEnd && this.prevFlowSignalsEnd < this.flowSignals.length) {
                    blockFlowSignals[idx][0] = new short[this.prevFlowSignalsEnd - this.prevFlowSignalsStart + 1];
                    if (this.readNegativeStrandFlag) {
                        for (j = this.prevFlowSignalsEnd; this.prevFlowSignalsStart <= j; j--) {
                            blockFlowSignals[idx][0][this.prevFlowSignalsEnd - j] = this.flowSignals[j];
                        }
                    } else {
                        for (j = this.prevFlowSignalsStart; j <= this.prevFlowSignalsEnd; j++) {
                            blockFlowSignals[idx][0][j - this.prevFlowSignalsStart] = this.flowSignals[j];
                        }
                    }
                } else {
                    blockFlowSignals[idx][0] = null;
                }
                // current context
                blockFlowSignals[idx][1] = new short[1];
                blockFlowSignals[idx][1][0] = s;
                // next context
                if (0 <= nextFlowSignalsStart && nextFlowSignalsStart <= nextFlowSignalsEnd && nextFlowSignalsEnd < this.flowSignals.length) {
                    blockFlowSignals[idx][2] = new short[nextFlowSignalsEnd - nextFlowSignalsStart + 1];
                    if (this.readNegativeStrandFlag) {
                        for (j = nextFlowSignalsEnd; nextFlowSignalsStart <= j; j--) {
                            blockFlowSignals[idx][2][nextFlowSignalsEnd - j] = this.flowSignals[j];
                        }
                    } else {
                        for (j = nextFlowSignalsStart; j <= nextFlowSignalsEnd; j++) {
                            blockFlowSignals[idx][2][j - nextFlowSignalsStart] = this.flowSignals[j];
                        }
                    }
                } else {
                    blockFlowSignals[idx][2] = null;
                }
                // update for the next iteration
                this.prevFlowSignalsStart = nextFlowSignalsStart;
                this.prevFlowSignalsEnd = nextFlowSignalsEnd;
                i++; // next base
                idx++; // next base
            }
            return blockFlowSignals;
        }
    }

    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     *
     * @param cigarString
     * @param readBases
     * @param readBaseQualities
     */
    private void createAlignmentBlocks(String cigarString, byte[] readBases, byte[] readBaseQualities) {
        createAlignmentBlocks(cigarString, readBases, readBaseQualities, null, null, -1);
    }

    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     *
     * @param cigarString
     * @param readBases
     * @param readBaseQualities
     * @param flowSignals       from the FZ tag, null if not present
     * @param flowOrder         from the RG.FO header tag, null if not present
     * @param flowOrderStart
     */
    private void createAlignmentBlocks(String cigarString, byte[] readBases, byte[] readBaseQualities,
                                       short[] flowSignals, String flowOrder, int flowOrderStart) {

        boolean showSoftClipped = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED);

        int nInsertions = 0;
        int nBlocks = 0;

        List<CigarOperator> operators = new ArrayList();
        StringBuffer buffer = new StringBuffer(4);

        if (cigarString.equals("*")) {
            alignmentBlocks = new AlignmentBlock[1];
            alignmentBlocks[0] = new AlignmentBlock(getStart(), readBases, readBaseQualities, this);
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
        FlowOrderAndSignalsBlockHelper fBlockHelper = null;
        if (null != flowSignals) {
            if (0 < readBases.length) {
                fBlockHelper = new FlowOrderAndSignalsBlockHelper(flowSignals, flowOrder, flowOrderStart, readBases, this.readNegativeStrandFlag);
            }
        }
        prevOp = 0;
        for (CigarOperator op : operators) {
            try {

                if (op.operator == HARD_CLIP) {
                    continue;
                }
                if (operatorIsMatch(showSoftClipped, op.operator)) {

                    byte[] blockBases = new byte[op.nBases];
                    byte[] blockQualities = new byte[op.nBases];
                    short[][][] blockFlowSignals = null;
                    AlignmentBlock block = null;

                    //Default value
                    Arrays.fill(blockQualities, (byte) 126);

                    int nBasesAvailable = readBases.length - fromIdx;

                    // TODO -- represent missing sequence ("*") explicitly for efficiency.
                    if (readBases == null || readBases.length == 0) {
                        Arrays.fill(blockBases, (byte) '=');
                    } else if (nBasesAvailable < op.nBases) {
                        Arrays.fill(blockBases, (byte) '?');
                    } else {
                        System.arraycopy(readBases, fromIdx, blockBases, 0, op.nBases);
                    }

                    nBasesAvailable = readBaseQualities.length - fromIdx;
                    if (readBaseQualities == null || readBaseQualities.length == 0 || nBasesAvailable < op.nBases) {
                        Arrays.fill(blockQualities, (byte) 126);
                    } else {
                        System.arraycopy(readBaseQualities, fromIdx, blockQualities, 0, op.nBases);
                    }

                    if (null != fBlockHelper) {
                        blockFlowSignals = fBlockHelper.createBlockFlowSignals(readBases, fromIdx, op.nBases);
                        block = AlignmentBlock.getInstance(blockStart, blockBases, blockQualities, blockFlowSignals, this);
                    } else {
                        block = AlignmentBlock.getInstance(blockStart, blockBases, blockQualities, this);
                    }

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
                    AlignmentBlock block = null;
                    short[][][] blockFlowSignals = null;

                    // This gap is between blocks split by insertion.   It is a zero
                    // length gap but must be accounted for.
                    gapTypes[gapIdx++] = ZERO_GAP;

                    byte[] blockBases = new byte[op.nBases];
                    byte[] blockQualities = new byte[op.nBases];

                    if (readBases == null || readBases.length == 0) {
                        Arrays.fill(blockBases, (byte) '=');
                    } else {
                        System.arraycopy(readBases, fromIdx, blockBases, 0, op.nBases);
                    }

                    if (readBaseQualities == null || readBaseQualities.length == 0) {
                        Arrays.fill(blockQualities, (byte) 126);
                    } else {
                        System.arraycopy(readBaseQualities, fromIdx, blockQualities, 0, op.nBases);
                    }

                    if (null != fBlockHelper) {
                        blockFlowSignals = fBlockHelper.createBlockFlowSignals(readBases, fromIdx, op.nBases);
                        block = AlignmentBlock.getInstance(blockStart, blockBases, blockQualities, blockFlowSignals, this);
                    } else {
                        block = AlignmentBlock.getInstance(blockStart, blockBases, blockQualities, this);
                    }

                    insertions[insertionIdx++] = block;

                    fromIdx += op.nBases;

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

    private boolean operatorIsMatch(boolean showSoftClipped, char operator) {
        return operator == MATCH || operator == PERFECT_MATCH || operator == MISMATCH
                || (showSoftClipped && operator == SOFT_CLIP);
    }

    public boolean isFirstInPair() {
        return firstInPair;
    }

    public boolean isNegativeStrand() {
        return readNegativeStrandFlag;
    }

    public boolean isDuplicate() {
        return duplicateReadFlag;
    }

    public boolean isMapped() {
        return !readUnmappedFlag;
    }

    public int getReadLength() {
        return this.readLength;
    }

    public boolean isPaired() {
        return readPairedFlag;
    }

    public boolean isProperPair() {
        return properPairFlag;
    }

    public boolean isFirstOfPair() {
        return this.firstRead;
    }

    public boolean isSecondOfPair() {
        return this.secondRead;
    }

    /**
     * Note: This method is required by the interface, but never used.
     *
     * @return
     */
    public LocusScore copy() {
        return new SamAlignment(record);
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

    public String getReadSequence() {
        return readSequence;
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

    public String getLibrary() {
        return library;
    }

    public SAMRecord getRecord() {
        return this.record;
    }

    @Override
    public String toString() {
        return record.getSAMString();
    }

    @Override
    public char[] getGapTypes() {
        return gapTypes;
    }

    public Object getAttribute(String key) {
        // SAM alignment tag keys must be of length 2
        return key.length() == 2 ? record.getAttribute(key) : null;
    }

    public String getClipboardString(double location) {
        return getValueStringImpl(location, false);
    }

    /**
     * Return info string for popup text and to copy to clipboard
     *
     * @param position
     * @param windowFunction -- not relevant, ignored
     * @return
     */
    public String getValueString(double position, WindowFunction windowFunction) {
        return getValueStringImpl(position, true);
    }

    String getValueStringImpl(double position, boolean truncate) {

        StringBuffer buf = new StringBuffer(super.getValueString(position, null));

        if (isPaired()) {
            boolean sectionBreak = false;
            if (record.getFirstOfPairFlag()) {
                buf.append("<br>First in pair");
                sectionBreak = true;
            }
            if (record.getSecondOfPairFlag()) {
                buf.append("<br>Second in pair");
                sectionBreak = true;
            }
            if (record.getNotPrimaryAlignmentFlag()) {
                buf.append("<br>Alignment NOT primary");
                sectionBreak = true;
            }
            if (record.getReadFailsVendorQualityCheckFlag()) {
                buf.append("<br>FAILED Vendor Quality Check");
                sectionBreak = true;
            }
            if (sectionBreak) {
                buf.append("<br>-------------------");
            }
        }

        List<SAMRecord.SAMTagAndValue> attributes = record.getAttributes();
        if (attributes != null && !attributes.isEmpty()) {

            for (SAMRecord.SAMTagAndValue tag : attributes) {
                buf.append("<br>" + tag.tag + " = ");

                if (tag.value.getClass().isArray()) { // ignore array types
                    buf.append("[not shown]<br>");
                    continue;
                }

                // Break tag
                final String tagValue = tag.value.toString();
                final int maxLength = 70;
                if (tagValue.length() > maxLength && truncate) {
                    String[] tokens = tagValue.split("<br>");
                    for (String token : tokens) {
                        if (token.length() > maxLength) {
                            // Insert line breaks
                            String remainder = token;
                            while (remainder.length() > maxLength) {
                                String tmp = remainder.substring(0, maxLength);
                                int spaceIndex = tmp.lastIndexOf(' ');
                                int idx = spaceIndex > 30 ? spaceIndex : maxLength;
                                final String substring = remainder.substring(0, idx);
                                buf.append(substring);
                                buf.append("<br>");
                                remainder = remainder.substring(idx);
                            }
                            buf.append(remainder);
                            buf.append("<br>");

                        } else {
                            buf.append(token);
                            buf.append("<br>");
                        }
                    }
                } else {
                    buf.append(tagValue);
                }

            }
            buf.append("<br>-------------------");
        }

        if (mateSequence != null) {
            buf.append("<br>Unmapped mate sequence: " + mateSequence);
            buf.append("<br>-------------------");
        }
        return buf.toString();
    }

    @Override
    public String getMateSequence() {
        return mateSequence;
    }

    public void setMateSequence(String mateSequence) {
        this.mateSequence = mateSequence;
    }

    public String getPairOrientation() {
        return pairOrientation;
    }

    public boolean isVendorFailedRead() {
        return record.getReadFailsVendorQualityCheckFlag();
    }

    public Color getDefaultColor() {
        return defaultColor;
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

    /**
     * @param record the SAM record
     * @return start index in the flow signal as specified by the ZF tag, or -1 if not present
     */
    public int getFlowSignalsStart(SAMRecord record) {
        Object attribute = record.getAttribute("ZF"); // NB: from a TMAP optional tag
        int toRet = -1;
        if(attribute != null && attribute instanceof Integer){
            toRet = (Integer) attribute;
        }
        return toRet;
    }

    /**
     * @param record      the SAM record
     * @param flowOrder   the flow order corresponding to this read
     * @param keySequence sequence the key sequence corresponding to this read
     * @return the flow signals in 100x format (SFF), only if they exist (FZ tag),
     *         if the key sequence and flow order are found in the read group header tag
     *         (RG.KS and RG.FO).  Note: the array proceeds in the sequencing direction.
     */
    public short[] getFlowSignals(SAMRecord record, String flowOrder, String keySequence) {
        short[] r = null;
        int i;
        int startFlow, keySignalOverlap;
        char firstBase;

        if (null == flowOrder || null == keySequence) {
            return null;
        }

        startFlow = this.getFlowSignalsStart(record);
        if (startFlow < 0) {
            return null;
        }

        // get the # of bases that the first base in the read overlaps with the last base(s) in the key
        if (this.readNegativeStrandFlag) {
            firstBase = (char) NT2COMP[record.getReadBases()[record.getReadLength() - 1]];
        } else {
            firstBase = (char) record.getReadBases()[0];
        }
        keySignalOverlap = 0;
        for (i = keySequence.length() - 1; 0 <= i && keySequence.charAt(i) == firstBase; i--) {
            keySignalOverlap += 100;
        }

        Object attribute = record.getAttribute("FZ");
        if (null == attribute) {
            return null;
        } else if (attribute instanceof short[]) {
            short[] signals = (short[]) attribute;
            r = new short[signals.length - startFlow];
            for (i = startFlow; i < signals.length; i++) {
                r[i - startFlow] = signals[i];
            }
        } else if (attribute instanceof int[]) {
            int[] signals = (int[]) attribute;
            r = new short[signals.length - startFlow];
            System.arraycopy(signals, startFlow, r, 0, r.length);
        } else if (attribute instanceof byte[]) {
            byte[] signals = (byte[]) attribute;
            r = new short[signals.length - startFlow];
            for (i = startFlow; i < signals.length; i++) {
                r[i - startFlow] = signals[i];
            }
        } else {
            return null;
        }
        // Subtract the key's contribution to the first base
        if (0 < keySignalOverlap && 0 < r.length) {
            if (r[0] <= keySignalOverlap) {
                r[0] = 0;
            } else {
                r[0] -= keySignalOverlap;
            }
        }

        return r;
    }
}
