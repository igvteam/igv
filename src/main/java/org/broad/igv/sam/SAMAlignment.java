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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
    public static final char UNKNOWN = 0;
    public static final String REDUCE_READS_TAG = "RR";

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
    protected int alignmentStart;
    protected int alignmentEnd;


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

    public SAMAlignment() {
    }


    public String getChr() {
        return chr;
    }

    @Override
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

    abstract public String getReadName();

    abstract public int getMappingQuality();

    abstract public int getInferredInsertSize();

    abstract public String getCigarString();

    abstract public String getReadLengthString();

    abstract public String getReadSequence();

    public AlignmentBlock[] getAlignmentBlocks() {
        return alignmentBlocks;
    }

    public AlignmentBlockImpl[] getInsertions() {
        return insertions;
    }

    public abstract boolean isNegativeStrand();

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
     */
    protected void createAlignmentBlocks(String cigarString, byte[] readBases, byte[] readBaseQualities) {

        if (cigarString.equals("*")) {
            alignmentBlocks = new AlignmentBlockImpl[1];
            alignmentBlocks[0] = new AlignmentBlockImpl(getStart(), readBases, readBaseQualities, readBases.length, '*');
            return;
        }

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
            gaps = new ArrayList<Gap>();
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

                    AlignmentBlockImpl block = buildAlignmentBlock(
                            op.operator,
                            readBases,
                            readBaseQualities,
                            blockStart,
                            fromIdx,
                            op.nBases,
                            true);

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
                    AlignmentBlockImpl block = buildAlignmentBlock(op.operator, readBases, readBaseQualities,
                            blockStart, fromIdx, op.nBases, false);
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


    /**
     * Build a list of cigar operators from a cigarString.  Removes padding operators and concatenates consecutive
     * operators of the same type
     *
     * @param cigarString
     * @return
     */
    public static List<CigarOperator> buildOperators(String cigarString) {

        java.util.List<CigarOperator> operators = new ArrayList();
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


    private static AlignmentBlockImpl buildAlignmentBlock(char operator,
                                                          byte[] readBases,
                                                          byte[] readBaseQualities,
                                                          int blockStart,
                                                          int fromIdx,
                                                          int nBases,
                                                          boolean checkNBasesAvailable) {

        byte[] blockBases = null;
        byte[] blockQualities = null;
        if (readBases != null && readBases.length > 0) {
            int nBasesAvailable = nBases;
            if (checkNBasesAvailable) {
                nBasesAvailable = readBases.length - fromIdx;
            }
            blockBases = new byte[nBases];
            if (nBasesAvailable < nBases) {
                Arrays.fill(blockBases, (byte) '?');
            }
            System.arraycopy(readBases, fromIdx, blockBases, 0, nBases);
        }
        if (readBaseQualities != null && readBaseQualities.length > 0) {
            int nBasesAvailable = nBases;
            if (checkNBasesAvailable) {
                nBasesAvailable = readBaseQualities.length - fromIdx;
            }
            blockQualities = new byte[nBases];
            if (nBasesAvailable < nBases) {
                Arrays.fill(blockQualities, (byte) 126);
            }
            System.arraycopy(readBaseQualities, fromIdx, blockQualities, 0, nBases);
        }
        AlignmentBlockImpl block = new AlignmentBlockImpl(blockStart, blockBases, blockQualities, nBases, operator);
        return block;
    }


    public String getClipboardString(double location, int mouseX) {
        return getValueStringImpl(location, mouseX, false);
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return getValueStringImpl(position, mouseX, true);
    }

    private String getValueStringImpl(double position, int mouseX, boolean truncate) {

        int basePosition = (int) position;
        StringBuffer buf = new StringBuffer();

        if (getHaplotypeName() != null) {
            buf.append("Hap name: " + getHaplotypeName() + "<br>");
            buf.append("Dist: " + getHapDistance() + "<br>");
        }

        // First check insertions.  Position is zero based, block coords 1 based
        if (this.insertions != null) {
            for (AlignmentBlock block : this.insertions) {

                if (block.containsPixel(mouseX)) {

                    byte[] bases = block.getBases();
                    if (bases == null) {
                        buf.append("Insertion: " + block.getLength() + " bases");
                    } else {
                        if (bases.length < 50) {
                            buf.append("Insertion (" + bases.length + " bases): " + new String(bases));
                        } else {
                            int len = bases.length;
                            buf.append("Insertion (" + bases.length + " bases): " + new String(Arrays.copyOfRange(bases, 0, 25)) + "..." +
                                    new String(Arrays.copyOfRange(bases, len - 25, len)));
                        }
                    }
                    return buf.toString();
                }
            }
        }

        // Not over an insertion

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


        String cigarString = getCigarString();
        // Abbreviate long CIGAR strings.  Retain the start and end of the CIGAR, which show
        // clipping; trim the middle.
        int maxCigarStringLength = 1000;
        if (cigarString.length() > maxCigarStringLength) {
            // Match only full <length><operator> pairs at the beginning and end of the string.
            Matcher lMatcher = Pattern.compile("^(.{1," + Integer.toString(maxCigarStringLength / 2 - 1) + "}[A-Z])").matcher(cigarString);
            Matcher rMatcher = Pattern.compile("[A-Z](.{1," + Integer.toString(maxCigarStringLength / 2) + "})$").matcher(cigarString);
            cigarString = (lMatcher.find() ? lMatcher.group(1) : "") + "..." + (rMatcher.find() ? rMatcher.group(1) : "");
        }


        buf.append("----------------------" + "<br>");
        buf.append("Mapping = " + (isPrimary() ? (isSupplementary() ? "Supplementary" : "Primary") : "Secondary") +
                (isDuplicate() ? " Duplicate" : "") + (isVendorFailedRead() ? " Failed QC" : "") +
                " @ MAPQ " + Globals.DECIMAL_FORMAT.format(getMappingQuality()) + "<br>");
        buf.append("Reference span = " + getChr() + ":" + Globals.DECIMAL_FORMAT.format(getAlignmentStart() + 1) + "-" +
                Globals.DECIMAL_FORMAT.format(getAlignmentEnd()) + " (" + (isNegativeStrand() ? "-" : "+") + ")" +
                " = " + Globals.DECIMAL_FORMAT.format(getAlignmentEnd() - getAlignmentStart()) + "bp<br>");
        buf.append("Cigar = " + cigarString + "<br>");
        buf.append("Clipping = ");

        // Identify the number of hard and soft clipped bases.
        Matcher lclipMatcher = Pattern.compile("^(([0-9]+)H)?(([0-9]+)S)?").matcher(cigarString);
        Matcher rclipMatcher = Pattern.compile("(([0-9]+)S)?(([0-9]+)H)?$").matcher(cigarString);
        int lclipHard = 0, lclipSoft = 0, rclipHard = 0, rclipSoft = 0;
        if (lclipMatcher.find()) {
            lclipHard = lclipMatcher.group(2) == null ? 0 : Integer.parseInt(lclipMatcher.group(2), 10);
            lclipSoft = lclipMatcher.group(4) == null ? 0 : Integer.parseInt(lclipMatcher.group(4), 10);
        }
        if (rclipMatcher.find()) {
            rclipHard = rclipMatcher.group(4) == null ? 0 : Integer.parseInt(rclipMatcher.group(4), 10);
            rclipSoft = rclipMatcher.group(2) == null ? 0 : Integer.parseInt(rclipMatcher.group(2), 10);
        }

        if (lclipHard + lclipSoft + rclipHard + rclipSoft == 0) {
            buf.append("None");
        } else {
            if (lclipHard + lclipSoft > 0) {
                buf.append("Left");
                if (lclipHard > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(lclipHard) + " hard");
                }
                if (lclipSoft > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(lclipSoft) + " soft");
                }
            }
            if (rclipHard + rclipSoft > 0) {
                buf.append((lclipHard + lclipSoft > 0 ? "; " : "") + "Right");
                if (rclipHard > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(rclipHard) + " hard");
                }
                if (rclipSoft > 0) {
                    buf.append(" " + Globals.DECIMAL_FORMAT.format(rclipSoft) + " soft");
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

        Object suppAlignment = this.getAttribute("SA");
        if (suppAlignment != null) {
            buf.append("----------------------<br>");
            buf.append(getSupplAlignmentString(suppAlignment.toString()));
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

                if (base == 0 && this.getReadSequence().equals("=") && !block.isSoftClipped() && genome != null) {
                    base = genome.getReference(chr, basePosition);

                }

                byte quality = block.getQuality(offset);
                buf.append("Location = " + getChr() + ":" + Globals.DECIMAL_FORMAT.format(1 + (long) position) + "<br>");
                buf.append("Base = " + (char) base + " @ QV " + Globals.DECIMAL_FORMAT.format(quality) + "<br>");

                break;
            }
        }

        return buf.toString();
    }


    // chr21,26002386,-,11785S1115M,60,0;chr21,26001844,+,1115S111M1D41M1D394M11239S,60,4;

    private String getSupplAlignmentString(String sa) {

        StringBuffer buf = new StringBuffer();
        buf.append("SupplementaryAlignments");
        String[] records = Globals.semicolonPattern.split(sa);
        for (String rec : records) {
            try {
                SupplementaryAlignment a = new SupplementaryAlignment(rec);
                buf.append("<br>" + a.printString());
            } catch (Exception e) {
                buf.append("<br>* Invalid SA entry (not listed) *");
            }
        }
        return buf.toString();
    }

    protected abstract String getAttributeString(boolean truncate);

    public abstract boolean isFirstOfPair();

    public abstract boolean isSecondOfPair();

    public abstract boolean isDuplicate();

    public abstract boolean isMapped();

    public abstract boolean isPaired();

    public abstract boolean isProperPair();

    public abstract boolean isSupplementary();

    public abstract boolean isVendorFailedRead();

    public abstract boolean isPrimary();

    /**
     * @return the unclippedStart
     */
    abstract public int getAlignmentStart();

    abstract public int getAlignmentEnd();

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


    public abstract Object getAttribute(String key);

    public java.util.List<Gap> getGaps() {
        return gaps;
    }


    @Override
    public void finish() {

    }

    @Override
    public AlignmentBlock getInsertionAt(int position) {
        for (AlignmentBlock block : insertions) {
            if (block.getStart() == position) return block;
            if (block.getStart() > position) return null;  // Blocks increase lineraly
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

    public static class SupplementaryAlignment {

        public String chr;
        public int start;
        public char strand;
        public int mapQ;
        public int numMismatches;
        public int lenOnRef;


        public SupplementaryAlignment(String rec) {
            String[] tokens = Globals.commaPattern.split(rec);
            chr = tokens[0];
            start = Integer.parseInt(tokens[1]);
            strand = tokens[2].charAt(0);
            mapQ = Integer.parseInt(tokens[4]);
            numMismatches = Integer.parseInt(tokens[5]);
            lenOnRef = computeLengthOnReference(tokens[3]);
        }

        public String printString() {
            // chr6:43,143,415-43,149,942 (-) @ MAPQ 60 NM 763
            return chr + ":" + Globals.DECIMAL_FORMAT.format(start) + "-" + Globals.DECIMAL_FORMAT.format(start + lenOnRef)
                    + " (" + strand + ") = " + Globals.DECIMAL_FORMAT.format(lenOnRef) + "bp  @MAPQ " + mapQ + " NM" + numMismatches;
        }


        int computeLengthOnReference(String cigarString) {

            int len = 0;
            StringBuffer buf = new StringBuffer();

            for (char c : cigarString.toCharArray()) {

                if (c > 47 && c < 58) {
                    buf.append(c);
                } else {
                    switch (c) {
                        case 'N':
                        case 'D':
                        case 'M':
                        case '=':
                        case 'X':
                            len += Integer.parseInt(buf.toString());
                    }
                    buf.setLength(0);
                }

            }
            return len;
        }
    }

    public String getSynopsisString() {

        char st = isNegativeStrand() ? '-' : '+';
        Object nm = getAttribute("NM");
        String numMismatches = nm == null ? "?" : nm.toString();
        int lenOnRef = getAlignmentEnd() - getAlignmentStart();

        int[] clipping = getClipping(getCigarString());
        String clippingString = "";
        if (clipping[0] + clipping[1] + clipping[2] + clipping[3] > 0) {
            if (clipping[0] > 0) clippingString += clipping[0] + "H";
            if (clipping[1] > 0) clippingString += clipping[1] + "S";
            clippingString += " ... ";
            if (clipping[3] > 0) clippingString += clipping[3] + "S";
            if (clipping[2] > 0) clippingString += clipping[2] + "H";
        }

        return chr + ":" + Globals.DECIMAL_FORMAT.format(getAlignmentStart()) + "-" +
                Globals.DECIMAL_FORMAT.format(getAlignmentEnd())
                + " (" + st + ") = " + Globals.DECIMAL_FORMAT.format(lenOnRef) + "BP  @MAPQ=" + getMappingQuality() +
                " NM=" + numMismatches + " CLIPPING=" + clippingString;

    }


    public static int[] getClipping(String cigarString) {
        // Identify the number of hard and soft clipped bases.
        Matcher lclipMatcher = Pattern.compile("^(([0-9]+)H)?(([0-9]+)S)?").matcher(cigarString);
        Matcher rclipMatcher = Pattern.compile("(([0-9]+)S)?(([0-9]+)H)?$").matcher(cigarString);
        int lclipHard = 0, lclipSoft = 0, rclipHard = 0, rclipSoft = 0;
        if (lclipMatcher.find()) {
            lclipHard = lclipMatcher.group(2) == null ? 0 : Integer.parseInt(lclipMatcher.group(2), 10);
            lclipSoft = lclipMatcher.group(4) == null ? 0 : Integer.parseInt(lclipMatcher.group(4), 10);
        }
        if (rclipMatcher.find()) {
            rclipHard = rclipMatcher.group(4) == null ? 0 : Integer.parseInt(rclipMatcher.group(4), 10);
            rclipSoft = rclipMatcher.group(2) == null ? 0 : Integer.parseInt(rclipMatcher.group(2), 10);
        }
        return new int[]{lclipHard, lclipSoft, rclipHard, rclipSoft};
    }


    ///// EXPERIMENTAL

    String haplotypeName;

    @Override
    public void setHaplotypeName(String hap) {
        haplotypeName = hap;
    }

    @Override
    public String getHaplotypeName() {
        return haplotypeName;
    }

    int hapDistance;

    @Override
    public void setHapDistance(int dist) {
        this.hapDistance = dist;
    }

    @Override
    public int getHapDistance() {
        return hapDistance;
    }

}
