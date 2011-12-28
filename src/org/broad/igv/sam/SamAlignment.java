/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.*;
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
    private String readGroup;
    private String library;
    private String sample;

	private boolean firstInPair;

    /**
     * Constructs ...
     *
     * @param record
     */
    public SamAlignment(SAMRecord record) {

        this.record = record;


        String refName = record.getReferenceName();

        Genome genome = Globals.isHeadless() ? null : IGV.getInstance().getGenomeManager().getCurrentGenome();
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
        this.readLength = record.getReadLength();
        this.firstInPair = record.getFirstOfPairFlag();

        setMatePair(record, genome);
        setPairOrientation(record);
        setFirstReadStrand(record);
        createAlignmentBlocks(record.getCigarString(), record.getReadBases(), record.getBaseQualities());

        SAMFileHeader header = record.getHeader();
        if (header != null) {
            readGroup = (String) record.getAttribute("RG");
            if (readGroup != null) {
                SAMReadGroupRecord rgRec = header.getReadGroup(readGroup);
                if (rgRec != null) {
                    String platform = rgRec.getPlatform();
                    sample = rgRec.getSample();
                    library = rgRec.getLibrary();

                }
            }
        }

        Object colorTag = record.getAttribute("YC");
        if (colorTag != null) {
            try {
                defaultColor = ColorUtilities.stringToColor(colorTag.toString());
            } catch (Exception e) {
                log.error("Error interpreting color tag: " + colorTag, e);
                defaultColor = AlignmentRenderer.grey1;
            }
        }
    }

    // Default constructor to support unit tests
    SamAlignment() {
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
        prevOp = 0;
        for (CigarOperator op : operators) {
            try {

                if (op.operator == HARD_CLIP) {
                    continue;
                }
                if (operatorIsMatch(showSoftClipped, op.operator)) {

                    byte[] blockBases = new byte[op.nBases];
                    byte[] blockQualities = new byte[op.nBases];


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

                    AlignmentBlock block = new AlignmentBlock(blockStart, blockBases, blockQualities, this);
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

                    insertions[insertionIdx++] = new AlignmentBlock(blockStart, blockBases, blockQualities, this);


                    fromIdx += op.nBases;

                }
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
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

    @Override
    public String toString() {
        return record.getSAMString();
    }

    @Override
    public char[] getGapTypes() {
        return gapTypes;
    }

    public Object getAttribute(String key) {
        return record.getAttribute(key);
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
