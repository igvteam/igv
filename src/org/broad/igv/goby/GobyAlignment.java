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
import java.util.Arrays;

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
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
    private AlignmentBlock[] block = new AlignmentBlock[1];
    private AlignmentBlock[] insertionBlock;
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
        //  buildBlocks(entry);

        block[0] = new AlignmentBlock(entry.getPosition(), buildBases(), buildQualities());
        insertionBlock = buildInsertions();

    }

    /**
     * Construct the AlignmentBlocks corresponding to insertions found in this goby Entry.
     *
     * @return
     */
    private AlignmentBlock[] buildInsertions() {
        int insertionCount = 0;
        for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {


            if (var.getFrom().length() < var.getTo().length()) {
                insertionCount++;
            }
        }
        AlignmentBlock[] result = new AlignmentBlock[insertionCount];
        if (insertionCount == 0) {
            return result;
        }

        int index = 0;
        for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {


            final String to = var.getTo();
            if (var.getFrom().length() < to.length()) {
                final char[] insertion = var.getTo().toCharArray();
                final byte[] insertedBytes = new byte[insertion.length];
                for (int j = 0; j < insertion.length; j++) {
                    insertedBytes[j] = (byte) insertion[j];
                }
                result[index++] = new AlignmentBlock(entry.getPosition() + var.getPosition(),
                        insertedBytes, to.getBytes());
            }
        }
        return result;
    }

    private byte[] buildQualities() {
        byte[] result = new byte[entry.getTargetAlignedLength()];
        final int length = result.length;
        for (int i = 0; i < length; i++) {
            result[i] = 40;
        }
        for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {

            final ByteString toQuality = var.getToQuality();
            for (int j = 0; j < toQuality.size(); j++) {
                final int offset = var.getPosition() + j - 1;
                if (offset >= length) break;
                result[offset] = toQuality.byteAt(j);
            }
        }
        return result;
    }

    /**
     * Reconstruct the read sequence, using the convention that '=' denotes a match to the reference.
     * <p/>
     * Conventions for storing sequence variations in Goby alignments are described <a
     * href="https://docs.google.com/document/d/1AVjhU23Ijowblb2qTGIqBGL7Nja8G3u7Gl2I7HVrRFY/edit?hl=en&authkey=CJeWuM8O#">here</a>
     *
     * @return the sequence of the read, with >< inserted around bases where an insertion occured in the read.
     */
    private byte[] buildBases() {
        final int length = entry.getTargetAlignedLength();
        byte[] result = new byte[length];
        Arrays.fill(result, (byte) '=');

        for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {
            final String from = var.getFrom();
            if (from.length() > 0 && from.charAt(0) == '-') {
                // Insertion
                final int offset = var.getPosition() - 1;
                if (offset >= length) break;
                result[offset] = (byte) '>';
                final int offsetPlusOne = offset + 1;
                if (offsetPlusOne >= length) break;
                result[offsetPlusOne] = (byte) '<';

            } else {
                final String to = var.getTo();
                for (int j = 0; j < to.length(); j++) {
                    final int offset = var.getPosition() + j - 1;
                    if (offset >= length) break;
                    result[offset] = (byte) to.charAt(j);
                }
            }
        }
        return result;
    }

    public void buildBlocks(Alignments.AlignmentEntry entry) {
        ObjectArrayList<AlignmentBlock> blocks = new ObjectArrayList<AlignmentBlock>();
        ObjectArrayList<AlignmentBlock> insertionBlocks = new ObjectArrayList<AlignmentBlock>();

        int start = entry.getPosition();
        ByteArrayList bases = new ByteArrayList();
        ByteArrayList scores = new ByteArrayList();

        int offset = 0;
        if (entry.getSequenceVariationsCount() > 0) {
            int firstVarPosition = entry.getSequenceVariations(0).getPosition();
            for (int i = offset; i < firstVarPosition; ++i) {
                bases.add((byte) '=');
                scores.add((byte) 40);
            }
        }
        for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {


            final String to = var.getTo();
            final ByteString toQuality = var.getToQuality();
            if (to.length() > 0 && to.charAt(0) == '-') {
                // insertion

                AlignmentBlock block = makeBlock(start, bases, scores);
                blocks.add(block);

                start += bases.size();
                offset = 0;
                bases.clear();
                scores.clear();

                for (int i = 0; i < to.length(); ++i) {
                    bases.add((byte) to.charAt(i));
                    if (i < toQuality.size()) {
                        scores.add(toQuality.byteAt(i));
                    }
                }

                AlignmentBlock insertionBlock = makeBlock(start, bases, scores);

                insertionBlocks.add(insertionBlock);
                bases.clear();
                scores.clear();
                // start is not changed when to = '--'
            } else {
                // other bases get appended (includes mutation and read insertion:
                for (int i = 0; i < to.length(); ++i) {
                    bases.add((byte) to.charAt(i));
                    if (i < toQuality.size()) {
                        scores.add(toQuality.byteAt(i));
                    }
                }


            }

        }
        for (int i = offset; i < entry.getQueryAlignedLength() - entry.getQueryPosition(); ++i) {
            // append bit that match after last variation until end of alignment:
            bases.add((byte) '=');
            scores.add((byte) 40);
        }
        AlignmentBlock oneBlock = makeBlock(start, bases, scores);
        blocks.add(oneBlock);
        start += bases.size();
        block = blocks.toArray(new AlignmentBlock[blocks.size()]);

        insertionBlock = insertionBlocks.toArray(new AlignmentBlock[insertionBlocks.size()]);
    }

    private AlignmentBlock makeBlock(int start, ByteArrayList bases, ByteArrayList scores) {
        return new AlignmentBlock(start,
                bases.toByteArray(new byte[bases.size()]),
                scores.toByteArray(new byte[scores.size()]));
    }

    static WarningCounter refMismatch = new WarningCounter(10);

    /**
     * Transform the read index into a readname:
     *
     * @return
     */
    public String getReadName() {
        //LOG.info("getReadName");
        return Integer.toString(entry.getQueryIndex());
    }

    public String getReadSequence() {
        //LOG.info("getReadSequence");
        return "read-sequence";
    }

    /**
     * Get the reference id from the iterator, prepend "chr".
     */
    public String getChromosome() {
        //LOG.info("getChromosome");
        return "chr" + iterator.indexToReferenceId.getId(entry.getTargetIndex()).toString();


    }

    /**
     * Get the reference id from the iterator, prepend "chr".
     */
    public String getChromosome(int targetIndex) {
        //LOG.info("getChromosome");
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

    public void setStart(int start) {
        //    //LOG.info("setStart");
        throw new UnsupportedOperationException("setStart is not supported");
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setEnd(int end) {
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
