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

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.mods.BaseModificationUtils;
import org.broad.igv.sam.mods.BaseModificationSet;
import org.broad.igv.sam.smrt.SMRTKinetics;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ultima.FlowUtil;
import org.broad.igv.ultima.annotate.FlowBlockAnnotator;

import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author jrobinso
 */
public class SAMAlignment implements Alignment {

    public static final int MAX_CIGAR_STRING_LENGTH_TO_DISPLAY = 50;
    public static final Pattern RIGHT_CIGAR_PATTERN = Pattern.compile("[A-Z](.{1," + MAX_CIGAR_STRING_LENGTH_TO_DISPLAY / 2 + "})$");
    public static final Pattern LEFT_CIGAR_PATTERN = Pattern.compile("^(.{1," + (MAX_CIGAR_STRING_LENGTH_TO_DISPLAY / 2 - 1) + "}[A-Z])");
    private static final Logger log = LogManager.getLogger(SAMAlignment.class);
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
    public static final char UNKNOWN = 0;
    public static final String REDUCE_READS_TAG = "RR";

    private static final int READ_PAIRED_FLAG = 0x1;
    private static final int PROPER_PAIR_FLAG = 0x2;
    private static final int READ_UNMAPPED_FLAG = 0x4;
    private static final int MATE_UNMAPPED_FLAG = 0x8;
    private static final int READ_STRAND_FLAG = 0x10;
    protected static final int MATE_STRAND_FLAG = 0x20;
    private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;
    private static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    private static final int DUPLICATE_READ_FLAG = 0x400;
    private static final int SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;

    private SAMReadGroupRecord readGroupRecord;

    private int flags;

    final private static FlowBlockAnnotator flowBlockAnnotator = new FlowBlockAnnotator();

    /**
     * Picard object upon which this SAMAlignment is based
     */
    private SAMRecord record;

    String chr;
    protected int start;  // <= Might differ from alignment start if soft clipping is considered
    protected int end;    // ditto
    protected Color ycColor = null;

    ReadMate mate;
    public AlignmentBlockImpl[] alignmentBlocks;
    public AlignmentBlockImpl[] insertions;
    List<Gap> gaps;
    char[] gapTypes;
    protected String mateSequence = null;
    protected String pairOrientation = "";
    private Strand firstOfPairStrand;
    private Strand secondOfPairStrand;

    /**
     * Map of position -> base modification
     */
    private List<BaseModificationSet> baseModificationSets;
    private SMRTKinetics smrtKinetics;

    private enum CacheKey {CLIPPING_COUNTS, SA_GROUP}

    ;

    /**
     * Use this to cache an expensive to compute value on this record and retrieve it as necessary.
     */
    @SuppressWarnings("unchecked")
    private <T> T getCachedOrCompute(CacheKey key, Supplier<T> supplier) {
        final Object value = record.getTransientAttribute(key);
        if (value == null) {
            final T newValue = supplier.get();
            record.setTransientAttribute(key, newValue);
            return newValue;
        } else {
            return (T) value;
        }
    }

    public SAMAlignment(SAMRecord record) {

        this.record = record;
        this.flags = record.getFlags();

        String refName = record.getReferenceName();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.chr = genome == null ? refName : genome.getCanonicalChrName(refName);

        // SAMRecord is 1 based inclusive.  IGV is 0 based exclusive.
        this.end = record.getAlignmentEnd();   // might be modified later for soft clipping
        this.start = record.getAlignmentStart() - 1;   // might be modified later for soft clipping

        if (record.getReadPairedFlag()) {
            String mateReferenceName = record.getMateReferenceName();
            String mateChr = genome == null ? mateReferenceName : genome.getCanonicalChrName(mateReferenceName);
            this.setMate(new ReadMate(mateChr,
                    record.getMateAlignmentStart() - 1,
                    record.getMateNegativeStrandFlag(),
                    record.getMateUnmappedFlag()));
        }

        SAMFileHeader header = record.getHeader();
        if (header != null) {
            String readGroup = (String) record.getAttribute("RG");
            if (readGroup != null) {
                this.readGroupRecord = header.getReadGroup(readGroup);

            }
        }

        Object colorTag = record.getAttribute("YC");
        if (colorTag != null) {
            try {
                ycColor = ColorUtilities.stringToColor(colorTag.toString(), null);
            } catch (Exception e) {
                log.error("Error interpreting color tag: " + colorTag, e);
            }
        }

        setPairOrientation();
        setPairStrands();
        createAlignmentBlocks();
    }

    public SAMRecord getRecord() {
        return this.record;
    }

    public String getChr() {
        return chr;
    }

    public String getContig() {
        return chr;
    }

    public String getDescription() {
        return getReadName();
    }

    public ReadMate getMate() {
        return mate;
    }

    public Color getYcColor() {
        return ycColor;
    }

    public Object getAttribute(String key) {
        // SAM alignment tag keys must be of length 2
        if (key == null) {
            return null;
        } else {
            return key.length() == 2 ? record.getAttribute(key) :
                    (key.equals("TEMPLATE_ORIENTATION") ? pairOrientation : null);
        }
    }

    public List<SAMRecord.SAMTagAndValue> getAttributes() {
        return record.getAttributes();
    }


    private Object getAttribute(SAMTag key) {
        return key == null ? null : record.getAttribute(key);
    }

    public boolean isFirstOfPair() {
        return isPaired() && (flags & FIRST_OF_PAIR_FLAG) != 0;
    }

    public boolean isSecondOfPair() {
        return isPaired() && (flags & SECOND_OF_PAIR_FLAG) != 0;
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

    public boolean isProperPair() {
        return ((flags & READ_PAIRED_FLAG) != 0) && ((flags & PROPER_PAIR_FLAG) != 0);
    }

    public boolean isNegativeStrand() {
        return (flags & READ_STRAND_FLAG) != 0;
    }

    public boolean isSupplementary() {
        return (flags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }

    public boolean isVendorFailedRead() {
        return (flags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    public boolean isPrimary() {
        return (flags & NOT_PRIMARY_ALIGNMENT_FLAG) == 0;
    }

    public String toString() {
        return record.getSAMString();
    }

    public String getReadName() {
        return record.getReadName();
    }

    public int getMappingQuality() {
        return record.getMappingQuality();
    }

    public int getInferredInsertSize() {
        return record.getInferredInsertSize();
    }

    public String getCigarString() {
        return record.getCigarString();
    }

    @Override
    public Cigar getCigar() {
        return record.getCigar();
    }

    @Override
    public ClippingCounts getClippingCounts() {
        return getCachedOrCompute(CacheKey.CLIPPING_COUNTS,
                () -> ClippingCounts.fromCigar(record.getCigar()));
    }

    public String getReadSequence() {
        return record.getReadString();
    }

    public int getAlignmentStart() {
        return record.getAlignmentStart() - 1;
    }

    public int getAlignmentEnd() {
        return record.getAlignmentEnd();
    }

    public String getReadLengthString() {

        Cigar cigar = getCigar();
        int readSequenceLength = cigar.getReadLength();
        int clippedLength = 0;
        CigarElement first = cigar.getFirstCigarElement();
        if (first.getOperator() == htsjdk.samtools.CigarOperator.H) {
            clippedLength += first.getLength();
        }
        CigarElement last = cigar.getLastCigarElement();
        if(last.getOperator() == htsjdk.samtools.CigarOperator.H) {
            clippedLength += last.getLength();
        }
        String readLengthString = Globals.DECIMAL_FORMAT.format(readSequenceLength + clippedLength) + " bp";
        if(clippedLength > 0) {
            readLengthString += " (" +  Globals.DECIMAL_FORMAT.format(readSequenceLength) + " sequence + " + Globals.DECIMAL_FORMAT.format(clippedLength) + " hard clipped)";
        }

        return readLengthString;
    }

    public String getSample() {
        return readGroupRecord == null ? null : readGroupRecord.getSample();
    }

    public String getReadGroup() {
        return readGroupRecord == null ? null : readGroupRecord.getId();
    }

    public String getLibrary() {
        return readGroupRecord == null ? null : readGroupRecord.getLibrary();
    }

    public AlignmentBlock[] getAlignmentBlocks() {
        return alignmentBlocks;
    }

    public AlignmentBlockImpl[] getInsertions() {
        return insertions;
    }

    public boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    public byte getBase(double position) {
        int basePosition = (int) position;
        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte base = block.getBase(offset);
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

    public List<BaseModificationSet> getBaseModificationSets() {

        if (baseModificationSets == null && record.hasAttribute("Mm") || record.hasAttribute("MM")) {

            Object mm = record.hasAttribute("Mm") ? record.getAttribute("Mm") : record.getAttribute("MM");
            byte[] ml = (byte[]) (record.hasAttribute("Ml") ? record.getAttribute("Ml") : record.getAttribute("ML"));

            // Minimal tag validation  -- 10X uses MM and/or ML for other purposes
            if (mm instanceof String && (mm.toString().length() > 0) && (ml == null || ml instanceof byte[])) {

                byte[] sequence = record.getReadBases();

                // Sequence length validation -- if MN tag is present use it, otherwise do a partial validation
                // Test sequence length vs mn if avaliable
                Integer mn = record.getIntegerAttribute("MN");
                if (mn == null) {
                    mn = sequence.length;
                }
                if (mn != null) {
                    if (mn != sequence.length) {
                        return null;
                    }
                } else if (!BaseModificationUtils.validateMMTag(this.getReadName(), mm.toString(), record.getReadBases(), isNegativeStrand())) {  //record.getCigarString().indexOf("H") > 0 &&
                    return null;
                }

                if (mm.toString().length() == 0) { // TODO -- more extensive validation?
                    baseModificationSets = Collections.EMPTY_LIST;
                } else {
                    baseModificationSets = BaseModificationUtils.getBaseModificationSets((String) mm, ml, sequence, isNegativeStrand());
                }
            }
        }
        return baseModificationSets;
    }

    public SMRTKinetics getSmrtKinetics() {
        if (smrtKinetics == null) {
            smrtKinetics = new SMRTKinetics(this);
        }
        return smrtKinetics;
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
                if (mate != null && mate.isMapped() && isProperPair()) {
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

    /**
     * Return the number of bases hard-clipped from the leading-edge (leftmost view in IGV) of the alignment
     */
    public int getLeadingHardClipLength() {
        int clipLength = 0;

        String cigarString = record.getCigarString();
        if (!cigarString.equals("*")) {
            java.util.List<CigarOperator> operators = buildOperators(cigarString);
            for (CigarOperator operator : operators) {
                if (operator.operator == HARD_CLIP) {
                    clipLength += operator.nBases;
                } else {
                    break;
                }
            }
        }
        return clipLength;
    }

    /**
     * Create the alignment blocks from the read bases and alignment information in the CIGAR
     * string.  The CIGAR string encodes insertions, deletions, skipped regions, and padding.
     */
    private void createAlignmentBlocks() {

        String cigarString = record.getCigarString();
        byte[] readBases = record.getReadBases();
        byte[] readBaseQualities = record.getBaseQualities();

        if (cigarString.equals("*")) {
            alignmentBlocks = new AlignmentBlockImpl[1];
            alignmentBlocks[0] = new AlignmentBlockImpl(getStart(), readBases, readBaseQualities, 0, readBases.length, '*');
        } else {
            // Create list of cigar operators
            java.util.List<CigarOperator> operators = buildOperators(cigarString);

            boolean showSoftClipped = PreferencesManager.getPreferences().getAsBoolean(Constants.SAM_SHOW_SOFT_CLIPPED);

            int nInsertions = 0;
            int nBlocks = 0;
            boolean firstOperator = true;
            int softClippedBaseCount = 0;
            int nGaps = 0;
            int nRealGaps = 0;
            char prevOp = 0;
            for (CigarOperator operator : operators) {

                char op = operator.operator;
                if (op == HARD_CLIP) {
                    continue;  // Just skip hardclips
                }
                int nBases = operator.nBases;
                if (operatorIsMatch(showSoftClipped, op)) {
                    nBlocks++;
                    if (operatorIsMatch(showSoftClipped, prevOp)) {
                        nGaps++;
                    }
                } else if (op == DELETION || op == SKIPPED_REGION) {
                    nGaps++;
                    nRealGaps++;
                } else if (op == INSERTION) {
                    nInsertions++;
                    nGaps++; // "virtual" gap, account for artificial block split @ insertion
                }

                if (firstOperator && op == SOFT_CLIP) {
                    softClippedBaseCount += nBases;
                }

                if (op != SOFT_CLIP) {
                    firstOperator = false;
                }

                prevOp = op;
            }


            alignmentBlocks = new AlignmentBlockImpl[nBlocks];
            insertions = new AlignmentBlockImpl[nInsertions];
            if (nGaps > 0) {
                gapTypes = new char[nGaps];
            }
            if (nRealGaps > 0) {
                gaps = new ArrayList<>();
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
            int padding = 0;
            prevOp = 0;
            for (int i = 0; i < operators.size(); i++) {
                CigarOperator op = operators.get(i);
                try {

                    if (op.operator == HARD_CLIP) {
                        continue;
                    }
                    if (operatorIsMatch(showSoftClipped, op.operator)) {

                        AlignmentBlockImpl block = AlignmentUtils.buildAlignmentBlock(
                                op.operator,
                                readBases,
                                readBaseQualities,
                                blockStart,
                                fromIdx,
                                op.nBases);

                        if (op.operator == SOFT_CLIP) {
                            block.setSoftClipped(true);
                        }
                        alignmentBlocks[blockIdx++] = block;

                        fromIdx += op.nBases;
                        blockStart += op.nBases;

                        if (operatorIsMatch(showSoftClipped, prevOp)) {
                            gapTypes[gapIdx++] = ZERO_GAP;
                        }

                    } else if (op.operator == DELETION) {
                        gaps.add(new Gap(blockStart, op.nBases, op.operator));
                        blockStart += op.nBases;
                        gapTypes[gapIdx++] = op.operator;
                    } else if (op.operator == SKIPPED_REGION) {

                        // Need the "flanking" regions, i.e. size of blocks either side of splice
                        // NOTE -- WE'RE ASSUMING HERE THAT THE "N" REGION IS FLANKED BY ALIGNMENT BLOCKS
                        int flankingLeft = 0;
                        int flankingRight = 0;
                        if (i > 0) {
                            flankingLeft = operators.get(i - 1).nBases;
                        }
                        if (i < operators.size() - 1) {
                            flankingRight = operators.get(i + 1).nBases;
                        }
                        gaps.add(new SpliceGap(blockStart, op.nBases, op.operator, flankingLeft, flankingRight));
                        blockStart += op.nBases;
                        gapTypes[gapIdx++] = op.operator;
                    } else if (op.operator == INSERTION) {
                        // This gap is between blocks split by insertion.   It is a zero
                        // length gap but must be accounted for.
                        gapTypes[gapIdx++] = ZERO_GAP;
                        AlignmentBlockImpl block = AlignmentUtils.buildAlignmentBlock(op.operator, readBases, readBaseQualities,
                                blockStart, fromIdx, op.nBases);
                        block.setPadding(padding);
                        insertions[insertionIdx++] = block;
                        fromIdx += op.nBases;
                        padding = 0;
                    } else if (op.operator == PADDING) {
                        // Padding for insertion start, which should be the next operator
                        padding += op.nBases;
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
    }


    /**
     * Build a list of cigar operators from a cigarString.  Removes padding operators and concatenates consecutive
     * operators of the same type
     *
     * @param cigarString
     * @return
     */
    private static List<CigarOperator> buildOperators(String cigarString) {

        java.util.List<CigarOperator> operators = new ArrayList<>();
        StringBuilder buffer = new StringBuilder(4);

        // Create list of cigar operators
        CigarOperator prevOp = null;
        for (int i = 0; i < cigarString.length(); i++) {
            char next = cigarString.charAt(i);
            if (Character.isDigit(next)) {
                buffer.append(next);
            } else {
                char op = next;
                int nBases = Integer.parseInt(buffer.toString());
                buffer.setLength(0);

                if (prevOp != null && prevOp.operator == op) {
                    prevOp.nBases += nBases;
                } else {
                    prevOp = new CigarOperator(nBases, op);
                    operators.add(prevOp);
                }
            }
        }
        return operators;

    }


    @Override
    public String getClipboardString(double location, int mouseX) {

        StringBuilder buf = new StringBuilder();

        String popupTextString = getAlignmentValueString(location, mouseX, null);
        buf.append(popupTextString);

        buf.append("----------------------");
        final String readSequence = getReadSequence();
        final String origSequence = isNegativeStrand() ? SequenceUtil.reverseComplement(readSequence) : readSequence;
        buf.append("\n\nRead sequence = ");
        buf.append(readSequence);
        buf.append("\n\nOrig sequence = ");
        buf.append(origSequence);
        buf.append("\n");

        return buf.toString();
    }

    private Integer positionToReadIndex(double position) {
        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains((int) position)) {
                return (int) (position - block.getStart()) + block.getBasesOffset();
            }
        }
        return null;
    }

    @Override
    public String getAlignmentValueString(double position, int mouseX, AlignmentTrack.RenderOptions renderOptions) {

        boolean truncate = renderOptions != null;
        int basePosition = (int) position;
        StringBuffer buf = new StringBuffer();

        if (getClusterName() != null) {
            buf.append("Cluster name: " + getClusterName() + "<br>");
            buf.append("Dist: " + getClusterDistance() + "<br>");
        }

        boolean hideSmallIndels = renderOptions == null ? false : renderOptions.isHideSmallIndels();
        int smallIndelThreshold = renderOptions == null ? 0 : renderOptions.getSmallIndelThreshold();

        boolean atInsertion = false;
        boolean atBaseMod = false;

        // First check insertions.  Position is zero based, block coords 1 based
        if (this.insertions != null) {
            for (AlignmentBlock block : this.insertions) {
                if (block.containsPixel(mouseX)) {
                    if (hideSmallIndels && block.getBasesLength() < smallIndelThreshold) {
                        continue;
                    }
                    ByteSubarray bases = block.getBases();
                    if (bases == null) {
                        buf.append("Insertion: " + block.getLength() + " bases<br>");
                    } else {
                        if (bases.length < 50) {
                            buf.append("Insertion (" + bases.length + " bases): " + bases.getString() + "<br>");
                        } else {
                            int len = bases.length;
                            buf.append("Insertion (" + bases.length + " bases): " + new String(bases.copyOfRange(0, 25)) + "..." +
                                    new String(bases.copyOfRange(len - 25, len)) + "<br>");
                        }

                        // extended annotation?
                        if (flowBlockAnnotator.handlesBlocks(block))
                            flowBlockAnnotator.appendBlockQualityAnnotation(this, block, buf);
                    }
                    atInsertion = true;
                }
            }
        }

        // Check base modifications & kinetics
        if (renderOptions != null) {
            final AlignmentTrack.ColorOption colorOption = renderOptions.getColorOption();
            if (colorOption.isBaseMod()) {
                Integer readIndex = positionToReadIndex(position);
                if (readIndex != null) {
                    int p = readIndex;
                    if (baseModificationSets != null) {
                        String modString = "";
                        for (BaseModificationSet bmSet : baseModificationSets) {
                            if (bmSet.containsPosition(p)) {
                                if (modString.length() > 0) modString += "<br>";
                                modString += bmSet.valueString(p);
                            }
                        }
                        if (modString.length() > 0) {
                            buf.append(modString);
                            buf.append("<br");
                            atBaseMod = true;
                        }
                    }
                }
            } else if (colorOption.isSMRTKinetics()) {
                Integer readIndex = positionToReadIndex(position);
                SMRTKinetics sk = getSmrtKinetics();
                if (readIndex != null) {
                    if (colorOption == AlignmentTrack.ColorOption.SMRT_SUBREAD_IPD) {
                        short[] ipdVals = sk.getSmrtSubreadIpd();
                        if (ipdVals != null) {
                            return "Subread IPD: " + Short.toUnsignedInt(ipdVals[readIndex]) + " Frames";
                        }
                    } else if (colorOption == AlignmentTrack.ColorOption.SMRT_SUBREAD_PW) {
                        short[] pwVals = sk.getSmrtSubreadPw();
                        if (pwVals != null) {
                            return "Subread PW: " + Short.toUnsignedInt(pwVals[readIndex]) + " Frames";
                        }
                    } else if (colorOption == AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD || colorOption == AlignmentTrack.ColorOption.SMRT_CCS_REV_IPD) {
                        final boolean isForwardStrand = (colorOption == AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD);
                        short[] ipdVals = sk.getSmrtCcsIpd(isForwardStrand);
                        if (ipdVals != null) {
                            final String strand = (isForwardStrand ? "fwd" : "rev");
                            return "CCS " + strand + "-strand aligned IPD: " + Short.toUnsignedInt(ipdVals[readIndex]) + " Frames";
                        }
                    } else {
                        final boolean isForwardStrand = (colorOption == AlignmentTrack.ColorOption.SMRT_CCS_FWD_PW);
                        short[] pwVals = sk.getSmrtCcsPw(isForwardStrand);
                        if (pwVals != null) {
                            final String strand = (isForwardStrand ? "fwd" : "rev");
                            return "CCS " + strand + "-strand aligned PW: " + Short.toUnsignedInt(pwVals[readIndex]) + " Frames";
                        }
                    }
                }
            }
        }

        if (atInsertion || atBaseMod) {
            return buf.toString();
        }

        buf.append("Read name = " + getReadName() + "<br>");

        String sample = getSample();
        if (sample != null) {
            buf.append("Sample = " + sample + "<br>");
        }
        String library = getLibrary();
        if (library != null) {
            buf.append("Library = " + library + "<br>");
        }
        String readGroup = getReadGroup();
        if (readGroup != null) {
            buf.append("Read group = " + readGroup + "<br>");
        }
        buf.append("Read length = " + getReadLengthString() + "<br>");

        buf.append("Flags = " + record.getFlags() + "<br>");

        buf.append("----------------------" + "<br>");
        String cigarString = getCigarString();
        // Abbreviate long CIGAR strings.  Retain the start and end of the CIGAR, which show
        // clipping; trim the middle.
        if (cigarString.length() > MAX_CIGAR_STRING_LENGTH_TO_DISPLAY) {
            // Match only full <length><operator> pairs at the beginning and end of the string.
            Matcher lMatcher = LEFT_CIGAR_PATTERN.matcher(cigarString);
            Matcher rMatcher = RIGHT_CIGAR_PATTERN.matcher(cigarString);
            cigarString = (lMatcher.find() ? lMatcher.group(1) : "") + "..." + (rMatcher.find() ? rMatcher.group(1) : "");
        }

        buf.append("Mapping = " + (isPrimary() ? (isSupplementary() ? "Supplementary" : "Primary") : "Secondary") +
                (isDuplicate() ? " Duplicate" : "") + (isVendorFailedRead() ? " Failed QC" : "") +
                " @ MAPQ " + Globals.DECIMAL_FORMAT.format(getMappingQuality()) + "<br>");
        buf.append("Reference span = " + getChr() + ":" + Globals.DECIMAL_FORMAT.format(getAlignmentStart() + 1) + "-" +
                Globals.DECIMAL_FORMAT.format(getAlignmentEnd()) + " (" + (isNegativeStrand() ? "-" : "+") + ")" +
                " = " + Globals.DECIMAL_FORMAT.format(getAlignmentEnd() - getAlignmentStart()) + "bp<br>");
        buf.append("Cigar = " + cigarString + "<br>");
        buf.append("Clipping = ");

        // Identify the number of hard and soft clipped bases.
        ClippingCounts clipping = getClippingCounts();

        if (!clipping.isClipped()) {
            buf.append("None");
        } else {
            if (clipping.isLeftClipped()) {
                buf.append("Left");
                if (clipping.getLeftHard() > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(clipping.getLeftHard()) + " hard");
                }
                if (clipping.getLeftSoft() > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(clipping.getLeftSoft()) + " soft");
                }
            }
            if (clipping.isRightClipped()) {
                buf.append((clipping.isLeftClipped() ? "; " : "") + "Right");
                if (clipping.getRightHard() > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(clipping.getRightHard()) + " hard");
                }
                if (clipping.getRightSoft() > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(clipping.getRightSoft()) + " soft");
                }
            }
        }
        buf.append("<br>");


        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        if (this.isPaired()) {
            buf.append("----------------------<br>");
            buf.append("Mate is mapped = " + (getMate().isMapped() ? "yes" : "no") + "<br>");
            if (getMate().isMapped()) {
                buf.append("Mate start = " + getMate().positionString() + "<br>");
                //buf.append("Pair is proper = " + (getProperPairFlag() ? "yes" : "no") + "<br>");
                if (getChr().equals(getMate().getChr())) {
                    buf.append("Insert size = " + getInferredInsertSize() + "<br>");
                }
            }
            if (isFirstOfPair()) {
                buf.append("First in pair<br>");
            }
            if (isSecondOfPair()) {
                buf.append("Second in pair<br>");
            }
            if (getPairOrientation().length() > 0) {
                buf.append("Pair orientation = " + getPairOrientation() + "<br>");
            }
        }

        Object suppAlignment = this.getAttribute(SAMTag.SA);
        if (suppAlignment != null) {
            buf.append("----------------------<br>");
            buf.append(getSupplAlignmentString());
            buf.append("<br>");
        }


        String attributeString = getAttributeString(truncate);
        if (attributeString != null && attributeString.length() > 0) {
            buf.append("----------------------");
            buf.append(getAttributeString(truncate));
        }


        if (mateSequence != null) {
            buf.append("----------------------<br>");
            buf.append("Mate sequence: " + mateSequence);
        }


        // Specific base

        for (AlignmentBlock block : this.alignmentBlocks) {
            if (block.contains(basePosition)) {

                buf.append("<hr>");
                int offset = basePosition - block.getStart();
                byte base = block.getBase(offset);

                if (base == 0 && this.getReadSequence().equals("=") && !block.isSoftClip() && genome != null) {
                    base = genome.getReference(chr, basePosition);

                }

                byte quality = block.getQuality(offset);
                buf.append("Location = " + getChr() + ":" + Globals.DECIMAL_FORMAT.format(1 + (long) position) + "<br>");
                buf.append("Base = " + (char) base + " @ QV " + Globals.DECIMAL_FORMAT.format(quality));
                if (FlowUtil.isFlow(this) && flowBlockAnnotator.handlesBlocks(block))
                    flowBlockAnnotator.appendBlockAttrAnnotation(this, block, offset, buf);
                buf.append("<br>");
                break;
            }
        }

        return buf.toString();
    }

    private String getAttributeString(boolean truncate) {
        // List of tags to skip.  Some tags, like MD and SA, are both quite verbose and not easily
        // interpreted by a human reader.  It is best to just hide these tags.  The list of tags
        // to hide is set through the SAM_HIDDEN_TAGS preference.
        ArrayList<String> tagsToHide = new ArrayList<String>(),
                tagsHidden = new ArrayList<String>();

        String samHiddenTagsPref = PreferencesManager.getPreferences().get(Constants.SAM_HIDDEN_TAGS);
        for (String s : (samHiddenTagsPref == null ? "" : samHiddenTagsPref).split("[, ]")) {
            if (!s.equals("")) {
                tagsToHide.add(s);
            }
        }

        StringBuffer buf = new StringBuffer();
        SAMRecord record = getRecord();
        List<SAMRecord.SAMTagAndValue> attributes = record.getAttributes();
        if (attributes != null && !attributes.isEmpty()) {

            for (SAMRecord.SAMTagAndValue tag : attributes) {
                if (tagsToHide.contains(tag.tag)) {
                    tagsHidden.add(tag.tag);
                    continue;
                }
                buf.append("<br>" + tag.tag + " = ");

                if (tag.tag.equals("ML") || tag.tag.equals("Ml")) {
                    buf.append(this.getMlTagString(tag));
                    buf.append("<br>");
                    continue;
                } else if (tag.value.getClass().isArray()) { // ignore array types
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

            if (tagsHidden.size() > 0) {
                buf.append("<br>Hidden tags: " + String.join(", ", tagsHidden));
            }
        }
        return buf.toString();
    }

    private String getMlTagString(SAMRecord.SAMTagAndValue tag) {
        byte[] bytes = (byte[]) tag.value;
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < bytes.length; i++) {
            if (i > 0) buf.append(",");
            buf.append(Byte.toUnsignedInt(bytes[i]));
        }
        return buf.toString();
    }

    private String getSupplAlignmentString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Supplementary Alignments");
        final List<SupplementaryAlignment> supplementaryAlignments = getSupplementaryAlignments();
        final int insertionIndex = SupplementaryAlignment.getInsertionIndex(this, supplementaryAlignments);
        int i = 0;
        for (SupplementaryAlignment sa : supplementaryAlignments) {
            if (i == insertionIndex) { //Add this read into the list
                sb.append(getThisReadDescriptionForSAList());
            }
            i++;
            try {
                sb.append("<br>" + sa.printString());
            } catch (Exception e) {
                sb.append("<br>* Invalid SA entry (not listed) *");
            }
        }
        if (i == insertionIndex) {
            sb.append(getThisReadDescriptionForSAList());
        }
        return sb.toString();
    }

    private String getThisReadDescriptionForSAList() {
        return "<br>"
                + "<b>"
                + chr + ":" + Globals.DECIMAL_FORMAT.format(getAlignmentStart() + 1) + "-" + Globals.DECIMAL_FORMAT.format(getAlignmentEnd())
                + " (" + this.getReadStrand().toShortString() + ") = " + Globals.DECIMAL_FORMAT.format(getLengthOnReference()) + "bp  @MAPQ " + this.getMappingQuality() + " NM" + getAttribute(SAMTag.NM)
                + "</b>";
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


    public java.util.List<Gap> getGaps() {
        return gaps;
    }

    @Override
    public AlignmentBlock getInsertionAt(int position) {
        for (AlignmentBlock block : insertions) {
            if (block.getStart() == position) return block;
            if (block.getStart() > position) return null;  // Blocks increase linearly
        }
        return null;
    }


    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public String getPairOrientation() {
        return pairOrientation;
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

    protected void setPairOrientation() {

        /*
                if (record.getReadPairedFlag() &&
                !record.getReadUnmappedFlag() &&
                !record.getMateUnmappedFlag() &&
                record.getReferenceName().equals(record.getMateReferenceName())) {

         */

        if (isPaired() && isMapped() && mate != null && mate.isMapped() && getChr().equals(mate.getChr())) {   // && name === mate.name

            char s1 = isNegativeStrand() ? 'R' : 'F';
            char s2 = mate.isNegativeStrand() ? 'R' : 'F';
            char o1 = ' ';
            char o2 = ' ';
            if (isFirstOfPair()) {
                o1 = '1';
                o2 = '2';
            } else if (isSecondOfPair()) {
                o1 = '2';
                o2 = '1';
            }

            final char[] tmp = new char[4];
            int isize = getInferredInsertSize();
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

    public void setChr(String chr) {
        this.chr = chr;
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

    public String getSynopsisString() {

        char st = isNegativeStrand() ? '-' : '+';
        Object nm = getAttribute(SAMTag.NM);
        String numMismatches = nm == null ? "?" : nm.toString();
        int lenOnRef = getAlignmentEnd() - getAlignmentStart();

        ClippingCounts clipping = getClippingCounts();
        String clippingString = "";
        if (clipping.isClipped()) {
            if (clipping.getLeftHard() > 0) clippingString += clipping.getLeftHard() + "H";
            if (clipping.getLeftSoft() > 0) clippingString += clipping.getLeftSoft() + "S";
            clippingString += " ... ";
            if (clipping.getRightSoft() > 0) clippingString += clipping.getRightSoft() + "S";
            if (clipping.getRightHard() > 0) clippingString += clipping.getRightHard() + "H";
        }

        return chr + ":" + Globals.DECIMAL_FORMAT.format(getAlignmentStart()) + "-" +
                Globals.DECIMAL_FORMAT.format(getAlignmentEnd())
                + " (" + st + ") = " + Globals.DECIMAL_FORMAT.format(lenOnRef) + "BP  @MAPQ=" + getMappingQuality() +
                " NM=" + numMismatches + " CLIPPING=" + clippingString;

    }

    private static boolean operatorIsMatch(boolean showSoftClipped, char operator) {
        return operator == MATCH || operator == PERFECT_MATCH || operator == MISMATCH
                || (showSoftClipped && operator == SOFT_CLIP);
    }


    /// // EXPERIMENTAL

    String haplotypeName;

    @Override
    public void setClusterName(String hap) {
        haplotypeName = hap;
    }

    @Override
    public String getClusterName() {
        return haplotypeName;
    }

    int hapDistance;

    @Override
    public void setHapDistance(int dist) {
        this.hapDistance = dist;
    }

    @Override
    public int getClusterDistance() {
        return hapDistance;
    }

    public List<SupplementaryAlignment> getSupplementaryAlignments() {
        return getCachedOrCompute(CacheKey.SA_GROUP,
                () -> {
                    Object rawSAValue = this.getAttribute(SAMTag.SA);
                    if (rawSAValue == null) {
                        return null;
                    } else {
                        final List<SupplementaryAlignment> sas = Collections.unmodifiableList(SupplementaryAlignment.parseFromSATag(rawSAValue.toString()));
                        return sas;
                    }
                });
    }
}
