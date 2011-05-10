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

package org.broad.igv.goby;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.EntryFlagHelper;
import edu.cornell.med.icb.goby.util.WarningCounter;
import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.track.WindowFunction;
import com.google.protobuf.ByteString;

import java.awt.*;
import java.util.*;

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import org.broad.igv.ui.IGV;

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

    private final Alignments.AlignmentEntry entry;
    private final GobyAlignmentIterator iterator;
    protected AlignmentBlock[] block = new AlignmentBlock[1];
    protected AlignmentBlock[] insertionBlock;
    private Color defaultColor = new Color(200, 200, 200);


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
        Arrays.fill(readQual, (byte) 40);
        int j = 0;

        final int leftPadding = alignmentEntry.getQueryPosition();
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
                    final byte qual = toQuality.size() > 0 ? toQuality.byteAt(i) : 40;

                    bases.add((byte) toChar);
                    scores.add(qual);

                }
                addBlock(insertionBlocks, alignmentEntry.getPosition() + var.getPosition(), bases, scores);
                bases.clear();
                scores.clear();
            } else if (!to.contains("-")) {
                for (int i = 0; i < toLength; i++) {
                    final int offset = j + var.getPosition() + i - 1 + leftPadding;
                    if (offset < readBases.length) {
                        readBases[offset] = (byte) to.charAt(i);
                    }
                }
            }
        }


        int pos = start;
        int endAlignmentRefPosition = readLength - leftPadding + start;
        bases.clear();
        scores.clear();
        while (pos < endAlignmentRefPosition) {

            bases.add(readBases[pos - start + leftPadding]);
            scores.add(readQual[pos - start + leftPadding]);
            ++pos;
        }

        addBlock(blocks, start, bases, scores);
        blocks = introduceDeletions(blocks, entry);
        block = blocks.toArray(new AlignmentBlock[blocks.size()]);

        insertionBlock = insertionBlocks.toArray(new AlignmentBlock[insertionBlocks.size()]);
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

                    AlignmentBlock left = new AlignmentBlock(block.getStart(),
                            leftBases.toByteArray(new byte[leftBases.size()]),
                            leftScores.toByteArray(new byte[leftScores.size()]));

                    AlignmentBlock right = new AlignmentBlock(block.getStart() + leftBases.size()
                            + var.getFrom().length(),
                            rightBases.toByteArray(new byte[rightBases.size()]),
                            rightScores.toByteArray(new byte[rightScores.size()]));

                    blocks.remove(block);
                    newBlocks.add(left);
                    newBlocks.add(right);

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

        blocks.add(new AlignmentBlock(start,
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
    public String getChromosome() {

        return "chr" + iterator.indexToReferenceId.getId(entry.getTargetIndex()).toString();


    }

    /**
     * Get the reference id from the iterator, prepend "chr".
     *
     * @param targetIndex Returns the chromosome id
     */
    public String getChromosome(int targetIndex) {

        return "chr" + iterator.indexToReferenceId.getId(targetIndex).toString();


    }

    public String getChr() {
        //LOG.info("getChr");
        return getChromosome();
    }

    public int getAlignmentStart() {
        //     //LOG.info("getAlignmentStart");
        return entry.getPosition();
    }

    public boolean contains(double location) {
        //LOG.info("contains");
        return false;  //To change body of implemented methods use File | Settings | File Templates.
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
        return new char[0];
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

    public ReadMate getMate() {
        if (entry.hasPairAlignmentLink()) {
            Alignments.RelatedAlignmentEntry link = entry.getPairAlignmentLink();
            String mateChr = getChromosome(link.getTargetIndex());
            int mateStart = entry.getPosition();
            boolean mateNegativeStrand = EntryFlagHelper.isMateReverseStrand(entry);

            boolean isReadUnmappedFlag = EntryFlagHelper.isReadUnmapped(entry);
            ReadMate mate = new ReadMate(mateChr, mateStart, mateNegativeStrand, isReadUnmappedFlag);
            return mate;
        } else {

            return null;  //To change body of implemented methods use File | Settings | File Templates.
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
    /*
public byte getBase(double position) {
//LOG.info("getBase");
return 0;
}

public byte getPhred(double position) {
for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {
for (int i = 0; i < var.getTo().length(); i++) {
 if (var.getPosition() + i == position) {
     return var.getToQuality().byteAt(i);
 }
}
}
// Goby only stores quality scores for variations at this point, so if we have not found the
// position in stored sequence variations, we return zero.
return 0;
}         */

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

    public Strand getFragmentStrand(int read) {
        //LOG.info("getFragmentStrand");
        return null;
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
    public Color getDefaultColor() {
        return defaultColor;
    }

    public int getStart() {
        //      //LOG.info("getStart");
        return entry.getPosition();
    }

    public int getEnd() {
        //    //LOG.info("getEnd");
        return entry.getPosition() + entry.getTargetAlignedLength();
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

    public String getValueString(double position, WindowFunction windowFunction) {
        //  //LOG.info("getValueString");
        MutableString buffer = new MutableString();

        buffer.append(entry.toString());
        buffer.replace("\n", "<br>");

        if (this.isPaired()) {
            buffer.append("----------------------" + "<br>");
            buffer.append("Pair start = " + getMate().positionString() + "<br>");
            buffer.append("Pair is mapped = " + (getMate().isMapped() ? "yes" : "no") + "<br>");
            //buf.append("Pair is proper = " + (getProperPairFlag() ? "yes" : "no") + "<br>");
            if (getChr().equals(getMate().getChr())) {
                buffer.append("Insert size = " + getInferredInsertSize() + "<br>");
            }
            if (getPairOrientation().length() > 0) {
                buffer.append("Pair orientation = " + getPairOrientation() + "<br>");
            }
        }
        return buffer.toString();
    }

    public byte getBase(double position) {
        int basePosition = (int) position;
        for (AlignmentBlock block : getAlignmentBlocks()) {
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
        for (AlignmentBlock block : getAlignmentBlocks()) {
            if (block.contains(basePosition)) {
                int offset = basePosition - block.getStart();
                byte score = block.getQuality(offset);
                return score;
            }
        }
        return 0;
    }
}
