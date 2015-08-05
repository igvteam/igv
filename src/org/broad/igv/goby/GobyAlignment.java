/*
 * The MIT License (MIT)
 *  Copyright (c) 2007-2015 by Institute for Computational Biomedicine,
 *                                          Weill Medical College of Cornell University.
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


package org.broad.igv.goby;

import com.google.protobuf.ByteString;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.EntryFlagHelper;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.CharArrayList;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.sam.*;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.Arrays;
import java.util.Comparator;

/**
 * A Facade to a <a href="http://goby.campagnelab.org">Goby</a> alignment entry. The facade exposes
 * <a href="http://goby.campagnelab.org">Goby</a> alignment entries in the format expected by
 * IGV. Since <a href="http://goby.campagnelab.org">Goby</a> does not store read sequences,
 * we retrieve the reference sequence on the fly from IGV and transform it to produce the read bases.
 * <p/>
 * For further information about Goby, or to obtain sample alignment files, see http://goby.campagnelab.org
 *
 * @author Fabien Campagne
 *         Date: Jun 29, 2010
 *         Time: 12:07:52 PM
 */
public class GobyAlignment implements Alignment {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(GobyAlignment.class);

    protected final Alignments.AlignmentEntry entry;
    private final GobyAlignmentIterator iterator;
    protected AlignmentBlock[] block = new AlignmentBlock[1];
    protected AlignmentBlock[] insertionBlock;
    private CharArrayList gapTypes = null;
    private static final ReadMate unmappedMate = new ReadMate("*", -1, false, true);
    private Comparator<? super AlignmentBlock> blockComparator = new Comparator<AlignmentBlock>() {
        public int compare(AlignmentBlock alignmentBlock, AlignmentBlock alignmentBlock1) {
            return alignmentBlock.getStart() - alignmentBlock1.getStart();
        }
    };

    /**
     * Construct the facade for an iterator and entry.
     *
     * @param iterator Used to retrieve chromosome names from target indices.
     * @param entry    Alignement entry (from Goby protocol buffer Alignment.entries).
     */
    public GobyAlignment(final GobyAlignmentIterator iterator, final Alignments.AlignmentEntry entry) {
        this.iterator = iterator;
        this.entry = entry;
        buildBlocks(entry);
    }


    private boolean hasReadInsertion(String from) {
        return from.length() > 0 && from.charAt(0) == '-';
    }

    /**
     * Construct alignment blocks from the Goby alignment entry. This method uses the convention that '=' denotes a match to the reference.
     * <p/>
     * Conventions for storing sequence variations in Goby alignments are described
     * <a href="http://tinyurl.com/goby-sequence-variations">here</a>
     *
     * @param alignmentEntry The Goby alignment entry to use
     */
    public void buildBlocks(Alignments.AlignmentEntry alignmentEntry) {

        ObjectArrayList<AlignmentBlock> blocks = new ObjectArrayList<AlignmentBlock>();
        ObjectArrayList<AlignmentBlock> insertionBlocks = new ObjectArrayList<AlignmentBlock>();

        int start = alignmentEntry.getPosition();
        ByteArrayList bases = new ByteArrayList();
        ByteArrayList scores = new ByteArrayList();
        int readLength = alignmentEntry.getQueryLength();

        byte[] readBases = new byte[readLength];
        byte[] readQual = new byte[readLength];
        Arrays.fill(readBases, (byte) '=');
        if (alignmentEntry.hasReadQualityScores()) {
            readQual = alignmentEntry.getReadQualityScores().toByteArray();
        } else {
            Arrays.fill(readQual, (byte) 40);
        }
        int j = 0;
        int insertedBases = 0;
        int deletedBases = 0;
        final int leftPadding = alignmentEntry.getQueryPosition();
        boolean showSoftClipped = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED);
        if (showSoftClipped && entry.hasSoftClippedBasesLeft()) {
            int clipLength = entry.getSoftClippedBasesLeft().length();

            addSoftClipBlock(blocks, Math.max(0,
                            entry.getPosition() - clipLength),
                    entry.getSoftClippedBasesLeft(), readQual,
                    entry.hasSoftClippedQualityLeft(),
                    entry.getSoftClippedQualityLeft().toByteArray(),
                    0);
        }
        for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
            final String from = var.getFrom();
            final int fromLength = from.length();
            final String to = var.getTo();
            final int toLength = from.length();
            final int sequenceVariationLength = Math.max(fromLength, toLength);
            final ByteString toQuality = var.getToQuality();

            if (hasReadInsertion(from)) {
                bases.clear();
                scores.clear();
                for (int i = 0; i < sequenceVariationLength; i++) {
                    final char toChar = i >= toLength ? '-' : to.charAt(i);
                    int size = toQuality.size();
                    final byte qual = size > 0 && i < size ? toQuality.byteAt(i) : 40;

                    bases.add((byte) toChar);
                    scores.add(qual);
                    deletedBases++;

                }
                addBlock(insertionBlocks, alignmentEntry.getPosition() + var.getPosition(), bases, scores);
                bases.clear();
                scores.clear();
            } else if (!to.contains("-")) {
                for (int i = 0; i < toLength; i++) {
                    final int offset = j + var.getPosition() + i - 1 + leftPadding - insertedBases;
                    if (offset > 0 && offset < readBases.length) {
                        readBases[offset] = (byte) to.charAt(i);
                        if (i < toQuality.size()) {
                            readQual[offset] = toQuality.byteAt(i);
                        }
                    }

                }
            } else {
                // has read deletion:
                insertedBases++;
            }
        }


        int pos = start;
        int matchLength = alignmentEntry.getQueryAlignedLength() - deletedBases;
        int endAlignmentRefPosition = matchLength + start;
        bases.clear();
        scores.clear();
        int maxIndex = Math.min(readBases.length, readQual.length);
        while (pos < endAlignmentRefPosition) {
            final int index = pos - start + leftPadding;
            if (index < maxIndex) {
                bases.add(readBases[index]);
                scores.add(readQual[index]);
            } else {
                break;
            }
            ++pos;
        }

        addBlock(blocks, start, bases, scores);
        blocks = introduceDeletions(blocks, entry);
        if (showSoftClipped && entry.hasSoftClippedBasesRight()) {

            int targetAlignedLength = entry.getTargetAlignedLength();
            addSoftClipBlock(blocks, entry.getPosition() + targetAlignedLength,
                    entry.getSoftClippedBasesRight(),
                    readQual,
                    entry.hasSoftClippedQualityRight(),
                    entry.getSoftClippedQualityRight().toByteArray(),
                    entry.getQueryAlignedLength() + entry.getSoftClippedBasesLeft().length());
        }
        block = blocks.toArray(new AlignmentBlock[blocks.size()]);
        Arrays.sort(block, blockComparator);
        insertionBlock = insertionBlocks.toArray(new AlignmentBlock[insertionBlocks.size()]);
        Arrays.sort(insertionBlock, blockComparator);
        ObjectArrayList<GobyAlignment> list = null;

        if (alignmentEntry.hasSplicedForwardAlignmentLink() || alignmentEntry.hasSplicedBackwardAlignmentLink()) {
            // if has a forward link, store a reference to this alignment in the reader (which represents the window scope)
            list = iterator.cacheSpliceComponent(this);
            if (list.size() > 1 && spliceListIsValid(list)) {

                final GobyAlignment spliceHeadAlignment = list.get(0);

                ObjectArrayList<AlignmentBlock> splicedBlocks = new ObjectArrayList<AlignmentBlock>();
                splicedBlocks.addAll(ObjectArrayList.wrap(spliceHeadAlignment.block));
                splicedBlocks.addAll(blocks);
                spliceHeadAlignment.block = splicedBlocks.toArray(new AlignmentBlock[splicedBlocks.size()]);

                ObjectArrayList<AlignmentBlock> splicedInsertionBlocks = new ObjectArrayList<AlignmentBlock>();
                splicedInsertionBlocks.addAll(ObjectArrayList.wrap(spliceHeadAlignment.insertionBlock));
                splicedInsertionBlocks.addAll(insertionBlocks);
                spliceHeadAlignment.insertionBlock = splicedInsertionBlocks.toArray(new AlignmentBlock[splicedInsertionBlocks.size()]);

                if (spliceHeadAlignment.gapTypes == null) {
                    spliceHeadAlignment.gapTypes = new CharArrayList(10);
                }
                spliceHeadAlignment.gapTypes.add(SAMAlignment.SKIPPED_REGION);

                // Since the previous alignment carries this information, we clear up block and insertionBlock
                // in this alignment, but keep any softClips:
                this.block = keepSoftClips(block);
                this.insertionBlock = new AlignmentBlock[0];
            }
        }

        block = removeNulls(block);

    }

    private AlignmentBlock[] removeNulls(AlignmentBlock[] block) {
        int nullCount = 0;
        for (int i = 0; i < block.length; i++) {
            AlignmentBlock alignmentBlock = block[i];
            if (alignmentBlock == null) {
                nullCount++;
            }
        }
        if (nullCount == 0) {
            // nothing to filter
            return block;
        } else {
            int newLength = block.length - nullCount;
            AlignmentBlock[] result = new AlignmentBlock[newLength];
            int j = 0;
            for (int i = 0; i < result.length; i++) {
                result[i] = block[j++];
            }
            return result;
        }
    }

    private AlignmentBlock[] keepSoftClips(AlignmentBlock[] blocks) {
        int numSoftCLippedBlocks = 0;
        for (AlignmentBlock block : blocks) {
            if (block.isSoftClipped()) numSoftCLippedBlocks++;
        }
        AlignmentBlock[] tmp = new AlignmentBlock[numSoftCLippedBlocks];
        int j = 0;
        for (int i = 0; i < numSoftCLippedBlocks; i++) {
            AlignmentBlock block = blocks[j++];
            if (block.isSoftClipped()) {
                tmp[i] = block;
            }
        }
        return tmp;
    }

    private void addSoftClipBlock(ObjectArrayList<AlignmentBlock> blocks, int position, String softClippedBasesLeft,
                                  byte[] readQualScores, boolean hasSoftClippedQuality,
                                  byte[] softClippedQuality, int j) {
        final int length = softClippedBasesLeft.length();
        byte[] bases = new byte[length];
        byte[] scores = new byte[length];

        for (int i = 0; i < length; i++) {
            bases[i] = (byte) softClippedBasesLeft.charAt(i);
            scores[i] = hasSoftClippedQuality ? softClippedQuality[i] : readQualScores[j++];
        }
        final AlignmentBlock alignmentBlock = new AlignmentBlock(getChr(), position,
                bases,
                scores);
        alignmentBlock.setSoftClipped(true);
        blocks.add(alignmentBlock);

    }

    /**
     * Verify that the list has an appropriate unbroken chain of back links.
     *
     * @param list the list of splices to validate
     * @return true if the list has an unbroken chain of back links
     */
    boolean spliceListIsValid(final ObjectArrayList<GobyAlignment> list) {

        if (list != null && list.size() > 1) {
            Alignments.AlignmentEntry prevEntry = list.get(0).entry;
            for (int i = 1; i < list.size(); i++) {
                Alignments.AlignmentEntry currentEntry = list.get(i).entry;
                if (!currentEntry.hasSplicedBackwardAlignmentLink()) return false;
                else {
                    Alignments.RelatedAlignmentEntry currentBackwardLink = currentEntry.getSplicedBackwardAlignmentLink();

                    if ((prevEntry.getQueryIndex() != currentEntry.getQueryIndex()) ||

                            (prevEntry.getFragmentIndex() != currentBackwardLink.getFragmentIndex()) ||
                            (prevEntry.getPosition() != currentBackwardLink.getPosition()) ||
                            (prevEntry.getTargetIndex() != currentBackwardLink.getTargetIndex())) {
                        return false;
                    }
                }
                prevEntry = currentEntry;
            }
        }
        return true;
    }


    /**
     * This method splits blocks whose boundaries contain a read deletion.
     *
     * @param blocks
     * @param alignmentEntry
     * @return
     */
    private ObjectArrayList<AlignmentBlock> introduceDeletions(ObjectArrayList<AlignmentBlock> blocks,
                                                               Alignments.AlignmentEntry alignmentEntry) {

        ObjectArrayList<AlignmentBlock> newBlocks = new ObjectArrayList<AlignmentBlock>();

        for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {


            for (AlignmentBlock block : blocks) {
                if (!block.isSoftClipped()) {

                    final int vrPos = var.getPosition() + entry.getPosition();
                    if (hasReadDeletion(var) && vrPos >= block.getStart() && vrPos <= block.getEnd()) {

                        ByteList leftBases = new ByteArrayList(block.getBases());
                        ByteList leftScores = new ByteArrayList(block.getQualities());
                        ByteList rightBases = new ByteArrayList(block.getBases());
                        ByteList rightScores = new ByteArrayList(block.getQualities());
                        int deletionPosition = var.getPosition() - 1;
                        leftBases = leftBases.subList(0, deletionPosition);
                        rightBases = rightBases.subList(deletionPosition, rightBases.size());

                        leftScores = leftScores.subList(0, deletionPosition);
                        rightScores = rightScores.subList(deletionPosition, rightScores.size());

                        AlignmentBlock left = new AlignmentBlock(getChr(), block.getStart(),
                                leftBases.toByteArray(new byte[leftBases.size()]),
                                leftScores.toByteArray(new byte[leftScores.size()]));

                        AlignmentBlock right = new AlignmentBlock(getChr(), block.getStart() + leftBases.size()
                                + var.getFrom().length(),
                                rightBases.toByteArray(new byte[rightBases.size()]),
                                rightScores.toByteArray(new byte[rightScores.size()]));

                        blocks.remove(block);
                        newBlocks.add(left);
                        newBlocks.add(right);

                    }
                }
            }
        }

        newBlocks.addAll(blocks);
        return newBlocks;
    }

    private boolean hasReadDeletion(Alignments.SequenceVariation var) {
        return (var.getTo().contains("-"));
    }


    private int addBlock(ObjectArrayList<AlignmentBlock> blocks, int start, ByteArrayList bases,
                         ByteArrayList scores) {

        blocks.add(
                new AlignmentBlock(getChr(), start,
                        bases.toByteArray(new byte[bases.size()]),
                        scores.toByteArray(new byte[scores.size()])));
        start += bases.size();
        bases.clear();
        scores.clear();
        return start;
    }


    /**
     * Transform the read index into a readname:
     *
     * @return
     */
    public String getReadName() {

        return Integer.toString(entry.getQueryIndex());
    }

    public String getReadSequence() {

        return "read-sequence";
    }

    /**
     * Get the reference id from the iterator, prepend "chr".
     */
    public String getChr() {
        return getChromosome(entry.getTargetIndex());
    }

    @Override
    public String getContig() {
        return getChr();
    }

    /**
     * Get the reference id from the iterator, prepend "chr".
     *
     * @param targetIndex Returns the chromosome id
     */
    public String getChromosome(int targetIndex) {
        return "chr" + iterator.getId(targetIndex).toString();
    }

    public int getAlignmentStart() {
        //     //LOG.info("getAlignmentStart");
        return entry.getPosition();
    }

    public boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }


    public AlignmentBlock[] getAlignmentBlocks() {
        //    //LOG.info("getAlignmentBlocks");

        return block;
    }

    public AlignmentBlock[] getInsertions() {
        //LOG.info("getInsertions");
        return insertionBlock;
    }


    public char[] getGapTypes() {
        //LOG.info("getGapTypes");
        if (gapTypes == null) {
            return new char[0];
        } else {
            return gapTypes.toArray();
        }
    }

    public String getCigarString() {
        //LOG.info("getCigarString");
        return null;
    }

    public int getInferredInsertSize() {
        if (entry.hasInsertSize()) {
            return entry.getInsertSize();
        } else return 0;
    }

    public int getMappingQuality() {
        if (entry.hasMappingQuality()) {
            return entry.getMappingQuality();
        } else {

            return 255;
        }
    }

    /**
     * Returns the mate for a paired-end read. Please note that this method will return an unmapped
     * mate for any single end read as well. Do check if the read is paired before calling getMate().
     *
     * @return The mate, or a constant unmapped mate (for single end reads, or paired end where the mate is not found).
     */
    public ReadMate getMate() {
        if (entry.hasPairAlignmentLink()) {
            Alignments.RelatedAlignmentEntry link = entry.getPairAlignmentLink();
            String mateChr = getChromosome(link.getTargetIndex());
            int mateStart = link.getPosition();
            boolean mateNegativeStrand = EntryFlagHelper.isMateReverseStrand(entry);

            boolean isReadUnmappedFlag = EntryFlagHelper.isReadUnmapped(entry);
            final ReadMate mate = new ReadMate(mateChr, mateStart, mateNegativeStrand, isReadUnmappedFlag);
            return mate;
        } else {
            return unmappedMate;
        }
    }

    public boolean isProperPair() {
        if (entry.hasPairFlags()) {

            return EntryFlagHelper.isProperlyPaired(entry);

        } else return false;
    }

    public boolean isMapped() {

        return true;
    }

    public boolean isPaired() {
        if (entry.hasPairFlags()) {

            return EntryFlagHelper.isPaired(entry);

        } else return false;
    }

    public boolean isNegativeStrand() {
        //     //LOG.info("isNegativeStrand");
        return entry.getMatchingReverseStrand();
    }

    public boolean isDuplicate() {
        //LOG.info("isDuplicate");
        return false;
    }

    public int getAlignmentEnd() {
        //LOG.info("getAlignmentEnd");
        return entry.getPosition() + entry.getTargetAlignedLength();
    }


    public String getSample() {
        //LOG.info("getSample");
        return null;
    }

    public String getReadGroup() {
        //LOG.info("getReadGroup");
        return null;
    }

    public String getLibrary() {
        //LOG.info("getReadGroup");
        return null;
    }

    public Object getAttribute(String key) {
        //LOG.info("getAttribute");
        return null;
    }

    /**
     * Return the strand of the read marked "first-in-pair" for a paired alignment. This method can return
     * Strand.NONE if the end marked first is unmapped.
     *
     * @return strand of first-of-pair
     */
    public Strand getFirstOfPairStrand() {
        if (isPaired()) {
            if (isFirstOfPair()) {
                return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
            } else {
                ReadMate mate = getMate();
                if (mate.isMapped() && isProperPair()) {
                    return mate.isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
                } else {
                    return Strand.NONE;
                }
            }

        } else {
            return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
        }
    }

    /**
     * Return the strand of the read marked "second-in-pair" for a paired alignment.  The strand is
     * undefined (Strand.NONE) for non-paired alignments
     *
     * @return strand of second-of-pair
     */
    public Strand getSecondOfPairStrand() {
        if (isPaired()) {
            if (isSecondOfPair()) {
                return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
            } else {
                ReadMate mate = getMate();
                if (mate.isMapped() && isProperPair()) {
                    return mate.isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
                } else {
                    return Strand.NONE;
                }
            }

        } else {
            // Undefined for non-paired alignments
            return Strand.NONE;
        }
    }

    public void setMateSequence(String sequence) {
        //LOG.info("setMateSequence");

    }

    public String getPairOrientation() {
        //LOG.info("getPairOrientation");
        String pairOrientation = "";
        if (EntryFlagHelper.isPaired(entry) &&

                !EntryFlagHelper.isMateUnmapped(entry) &&
                entry.getTargetIndex() == entry.getPairAlignmentLink().getTargetIndex()) {

            char s1 = EntryFlagHelper.isReadReverseStrand(entry) ? 'R' : 'F';
            char s2 = EntryFlagHelper.isMateReverseStrand(entry) ? 'R' : 'F';
            char o1 = ' ';
            char o2 = ' ';
            char[] tmp = new char[4];
            if (EntryFlagHelper.isFirstInPair(entry)) {
                o1 = '1';
                o2 = '2';
            } else if (EntryFlagHelper.isSecondInPair(entry)) {
                o1 = '2';
                o2 = '1';
            }
            if (getInferredInsertSize() > 0) {
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
        return pairOrientation;
    }

    public boolean isSmallInsert() {
        //LOG.info("isSmallInsert");
        return false;
    }


    /**
     * Return true if this read failed vendor quality checks
     */
    public boolean isVendorFailedRead() {
        return false;
    }

    /**
     * Return the default color with which to render this alignment
     *
     * @return
     */
    public Color getColor() {
        return null;
    }

    public int getStart() {
        //      //LOG.info("getStart");
        return entry.getPosition();
    }

    public int getEnd() {
        //    //LOG.info("getEnd");
        if (block.length == 0) return getStart();
        else {

            int last = block.length - 1;
            if (block[last] == null) {
                // System.out.println("STOP");
                return entry.getPosition() + entry.getTargetAlignedLength();
            }
            return block[last].getEnd();
        }
    }

    public void setStart(
            int start) {
        //    //LOG.info("setStart");
        throw new UnsupportedOperationException("setStart is not supported");
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setEnd(
            int end) {
        throw new UnsupportedOperationException("setEnd is not supported");
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public float getScore() {
        //LOG.info("getScore");
        return entry.getScore();
    }


    public LocusScore copy() {
        return this;
    }


    /**
     * This method is provide to provide as a hook for customizing the text that is copied to the clipboard.
     * The default behavior is to just copy the tooltip text.
     *
     * @param location
     * @return
     */
    public String getClipboardString(double location) {
        return getValueString(location, null);
    }


    public String getValueString(double position, WindowFunction windowFunction) {
        //  //LOG.info("getValueString");
        MutableString buffer = new MutableString();

        buffer.append(entry.toString());
        buffer.replace("\n", "<br>");

        if (this.isPaired()) {
            buffer.append("----------------------" + "<br>");
            buffer.append("Mate start = " + getMate().positionString() + "<br>");
            buffer.append("Mate is mapped = " + (getMate().isMapped() ? "yes" : "no") + "<br>");
            //buf.append("Pair is proper = " + (getProperPairFlag() ? "yes" : "no") + "<br>");
            if (getChr().equals(getMate().getChr())) {
                buffer.append("Insert size = " + getInferredInsertSize() + "<br>");
            }
            if (getPairOrientation().length() > 0) {
                buffer.append("Pair orientation = " + getPairOrientation() + "<br>");
            }
            if (isFirstOfPair()) {
                buffer.append("First of pair <br>");
            }
            if (isSecondOfPair()) {
                buffer.append("Second of pair <br>");
            }
        }
        return buffer.toString();
    }

    public byte getBase(double position) {
        int basePosition = (int) position;
        for (AlignmentBlock block : getAlignmentBlocks()) {
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
        for (AlignmentBlock block : getAlignmentBlocks()) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte score = block.getQuality(offset);
                return score;
            }
        }
        return 0;
    }


    /**
     * Return true if this alignment is marked "first in pair".  Added to suppor bisulfite sequencing mode.
     */
    public boolean isFirstOfPair() {
        return EntryFlagHelper.isFirstInPair(entry);
    }


    /**
     * Return true if this alignment is marked "second in pair".  Added to suppor bisulfite sequencing mode.
     */
    public boolean isSecondOfPair() {
        return EntryFlagHelper.isSecondInPair(entry);
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public void finish() {
    }

    @Override
    public boolean isPrimary() {
        return !EntryFlagHelper.isNotPrimaryAlignment(entry);
    }

    @Override
    public boolean isSupplementary() {
        // The SAM 0x0800 tag
        return false;
    }
}
