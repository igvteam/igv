/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.ColorUtilities;

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
    boolean negativeStrand;
    boolean readNegativeStrandFlag;
    boolean duplicateReadFlag;
    boolean readUnmappedFlag;
    boolean readPairedFlag;
    boolean properPairFlag;
    SAMRecord record;
    String cigarString;
    String readSequence;
    private boolean softClippedStart = false;
    private boolean softClippedEnd = false;
    private Strand firstReadStrand = Strand.NONE;
    private Strand secondReadStrand = Strand.NONE;
    boolean firstRead = false;
    boolean secondRead = false;
    private String mateSequence = null;
    private String pairOrientation = "";
    private Color defaultColor = AlignmentRenderer.grey1;


    /**
     * Constructs ...
     *
     * @param record
     */
    public SamAlignment(SAMRecord record) {

        this.record = record;


        String refName = record.getReferenceName();

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.chr = genome == null ? refName : genome.getChromosomeAlias(refName);


        // SAMRecord is 1 based inclusive.  IGV is 0 based exclusive.
        this.alignmentStart = record.getAlignmentStart() - 1;
        this.start = this.alignmentStart;   // might be modified later for soft clipping
        this.alignmentEnd = Math.max(alignmentStart, record.getAlignmentEnd());
        this.end = alignmentEnd;   // might be modified later for soft clipping
        this.negativeStrand = record.getReadNegativeStrandFlag();
        this.cigarString = record.getCigarString();
        this.setMappingQuality(record.getMappingQuality());
        this.readName = record.getReadName().trim();
        this.readNegativeStrandFlag = record.getReadNegativeStrandFlag();
        this.duplicateReadFlag = record.getDuplicateReadFlag();
        this.readUnmappedFlag = record.getReadUnmappedFlag();
        this.readPairedFlag = record.getReadPairedFlag();
        this.setInferredInsertSize(record.getInferredInsertSize());
        this.readSequence = record.getReadString();

        setMatePair(record, genome);
        setPairOrientation(record);
        setFirstReadStrand(record);
        createAlignmentBlocks(record.getCigarString(), record.getReadBases(), record.getBaseQualities());


        Object colorTag = record.getAttribute("YC");
        if(colorTag != null) {
            try {
                defaultColor = ColorUtilities.convertRGBStringToColor(colorTag.toString());
            } catch (Exception e) {
                log.error("Error interpreting color tag: " + colorTag, e);
                defaultColor = AlignmentRenderer.grey1;
            }
        }
    }

    private void setFirstReadStrand(SAMRecord record) {
        if (!record.getReadPairedFlag()) {
            firstReadStrand = (record.getReadNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE);
        } else if (record.getProperPairFlag()) {
            if (record.getFirstOfPairFlag()) {
                firstReadStrand = (record.getReadNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE);
                if (!record.getMateUnmappedFlag()) {
                    secondReadStrand = (record.getMateNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE);
                }
            } else {
                if (!record.getMateUnmappedFlag()) {
                    firstReadStrand = (record.getMateNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE);
                }
                secondReadStrand = (record.getReadNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE);
            }
        }
    }

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

    private char[] tmp = new char[4];

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
     * A copy constructor
     */
    private SamAlignment(SamAlignment alignment) {
        this.chr = alignment.chr;
        this.alignmentStart = alignment.alignmentStart;
        this.alignmentEnd = alignment.alignmentEnd;
        this.end = alignment.end;
        this.negativeStrand = alignment.negativeStrand;
        this.mate = alignment.mate;
        this.alignmentBlocks = alignment.alignmentBlocks;
        this.insertions = alignment.insertions;
        this.cigarString = alignment.cigarString;
        this.mappingQuality = alignment.mappingQuality;
        this.readName = alignment.readName;
        this.readNegativeStrandFlag = alignment.readNegativeStrandFlag;
        this.duplicateReadFlag = alignment.duplicateReadFlag;
        this.readUnmappedFlag = alignment.readUnmappedFlag;
        this.readPairedFlag = alignment.readPairedFlag;
        this.inferredInsertSize = alignment.inferredInsertSize;
        this.properPairFlag = alignment.properPairFlag;
    }
    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     *
     * @param cigarString
     * @param readBases
     * @param readBaseQualities
     */
    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     *
     * @param cigarString
     * @param readBases
     * @param readBaseQualities
     */

    void createAlignmentBlocks(String cigarString, byte[] readBases, byte[] readBaseQualities) {

        boolean showSoftClipped = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED); 


        int nInsertions = 0;
        int nBlocks = 0;

        List<CigarOperator> operators = new ArrayList();
        StringBuffer buffer = new StringBuffer(4);

        if (cigarString.equals("*")) {
            alignmentBlocks = new AlignmentBlock[1];
            alignmentBlocks[0] = new AlignmentBlock(getStart(), readBases, readBaseQualities);
            return;
        }

        boolean firstOperator = true;
        int softClippedBaseCount = 0;
        int nGaps = 0;
        char prevOp = 0;
        for (int i = 0; i < cigarString.length(); i++) {
            char next = cigarString.charAt(i);
            if (Character.isDigit(next)) {
                buffer.append(next);
            } else {
                int nBases = Integer.parseInt(buffer.toString());
                char op = next;
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
        if (nGaps > 0)

        {
            gapTypes = new char[nGaps];
        }

        // Adjust start and fill new arrays
        if (showSoftClipped) {
            start -= softClippedBaseCount;
        }
        int fromIdx = showSoftClipped ? 0 : softClippedBaseCount;
        int blockStart = start;

        int blockIdx = 0;
        int insertionIdx = 0;
        int gapIdx = 0;
        prevOp = 0;
        for (CigarOperator op : operators) {
            if (operatorIsMatch(showSoftClipped, op.operator)) {

                byte[] blockBases = new byte[op.nBases];
                byte[] blockQualities = new byte[op.nBases];
                //Default value
                Arrays.fill(blockQualities, (byte) 126);

                // TODO -- represent missing sequence ("*") explicitly rather.
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

                AlignmentBlock block = new AlignmentBlock(blockStart, blockBases, blockQualities);
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

                // This gap is between blocks split by insertion.   It is posA zero
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

                insertions[insertionIdx++] = new AlignmentBlock(blockStart, blockBases, blockQualities);

                fromIdx += op.nBases;

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

    // Default constructor to support unit tests

    SamAlignment() {
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

    public boolean isPaired() {
        return readPairedFlag;
    }

    public boolean isProperPair() {
        return properPairFlag;
    }


    /**
     * TODO
     *
     * @return
     */
    public LocusScore copy() {
        return new SamAlignment(this);
    }


    /**
     * @return the unclippedStart
     */
    public int getAlignmentStart() {
        return alignmentStart;
    }

    /**
     * @param unclippedStart the unclippedStart to set
     */
    public void setUnclippedStart(int unclippedStart) {
        this.alignmentStart = unclippedStart;
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
        this.alignmentStart = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
        this.alignmentEnd = end;
    }

    /**
     * @return the softClippedStart
     */
    public boolean isSoftClippedStart() {
        return softClippedStart;
    }

    /**
     * @return the softClippedEnd
     */
    public boolean isSoftClippedEnd() {
        return softClippedEnd;
    }

    public String getSample() {

        SAMFileHeader header = record.getHeader();
        if (header != null) {
            String rg = (String) record.getAttribute("RG");
            if (rg != null) {
                SAMReadGroupRecord rgRec = header.getReadGroup(rg);
                if (rgRec != null) {
                    return rgRec.getSample();
                }
            }
        }
        return null;
    }

    public String getReadGroup() {

        SAMFileHeader header = record.getHeader();
        if (header != null) {
            String readGroup = (String) record.getAttribute("RG");
            return (String) record.getAttribute("RG");
        }
        return null;
    }

    @Override
    public String toString() {
        return record.format();
    }

    @Override
    public char[] getGapTypes() {
        return gapTypes;
    }

    public Strand getFragmentStrand(int strand) {
        return strand == 1 ? firstReadStrand : secondReadStrand;
    }

    public Object getAttribute(String key) {
        return record.getAttribute(key);
    }

    /**
     * Return info string for popup text and to copy to clipboard
     *
     * @param position
     * @param windowFunction -- not relevant, ignored
     * @return
     */
    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer(super.getValueString(position, windowFunction));

        if (isPaired()) {
            if (record.getFirstOfPairFlag()) {
                buf.append("<br>First in pair");
            }
            if (record.getSecondOfPairFlag()) {
                buf.append("<br>Second in pair");
            }
            if (record.getNotPrimaryAlignmentFlag()) {
                buf.append("<br>Alignment NOT primary");
            }
            if (record.getReadFailsVendorQualityCheckFlag()) {
                buf.append("<br>FAILED Vendor Quality Check");
            }
            buf.append("<br>-------------------");
        }

        List<SAMRecord.SAMTagAndValue> attributes = record.getAttributes();
        if (attributes != null && !attributes.isEmpty()) {

            for (SAMRecord.SAMTagAndValue tag : attributes) {
                buf.append("<br>" + tag.tag + " = " + tag.value);
            }
            buf.append("<br>-------------------");
        }

        if (mateSequence != null) {
            buf.append("<br>Unmapped mate sequence: " + mateSequence);
            buf.append("<br>-------------------");
        }
        return buf.toString();

    }

    public boolean isFirstInPair() {
        return record.getFirstOfPairFlag();
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
