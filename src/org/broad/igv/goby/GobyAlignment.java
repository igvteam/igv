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
    private Color defaultColor = new Color(200, 200, 200);;

    /**
     * Construct the facade for an iterator and entry.
     *
     * @param iterator Used to retrieve chromosome names from target indices.
     * @param entry    Alignement entry (from Goby protocol buffer Alignment.entries).
     */
    public GobyAlignment(final GobyAlignmentIterator iterator, final Alignments.AlignmentEntry entry) {
        this.iterator = iterator;
        this.entry = entry;

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

    private byte[] buildBases() {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        String genomeId = genome.getId();

        String reference = iterator.getReference();
        String referenceAlias = genome.getChromosomeAlias(reference);
        int position = entry.getPosition();



        // adjust by one because Goby positions start at 1 while IGV starts at 0
        byte[] result = SequenceManager.readSequence(genomeId, referenceAlias, position, position + entry.getTargetAlignedLength());
        final int length = result.length;

        for (Alignments.SequenceVariation var : entry.getSequenceVariationsList()) {
            final String to = var.getTo();
            for (int j = 0; j < to.length(); j++) {

                final int offset = var.getPosition() + j - 1;
                if (offset >= length) break;
                result[offset] = (byte) to.charAt(j);
            }
        }
        return result;
    }

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
        return new char[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getCigarString() {
        //LOG.info("getCigarString");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getInferredInsertSize() {
        //LOG.info("getInferredInsertSize");
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getMappingQuality() {
        //LOG.info("getMappingQuality");
        return 255;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public ReadMate getMate() {
        //LOG.info("getMate");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isProperPair() {
        //LOG.info("isProperPair");
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isMapped() {
        //LOG.info("isMapped");
        return true;
    }

    public boolean isPaired() {
        //LOG.info("isPaired");
        return false;
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

    public byte getBase(double position) {
        //LOG.info("getBase");
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public byte getPhred(double position) {
        //LOG.info("getPhred");
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getSample() {
        //LOG.info("getSample");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getReadGroup() {
        //LOG.info("getReadGroup");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Object getAttribute(String key) {
        //LOG.info("getAttribute");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Strand getFragmentStrand(int read) {
        //LOG.info("getFragmentStrand");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setMateSequence(String sequence) {
        //LOG.info("setMateSequence");
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getPairOrientation() {
        //LOG.info("getPairOrientation");
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isSmallInsert() {
        //LOG.info("isSmallInsert");
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


    /**
     * Return true if this read failed vendor quality checks
      */    
    public boolean isVendorFailedRead() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
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
        return 1;
    }

    public void setConfidence(float confidence) {
        throw new UnsupportedOperationException("setEnd is not supported");
    }

    public float getConfidence() {
        //LOG.info("getConfidence");
        return entry.getScore();
    }

    public LocusScore copy() {
        return this;
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        //  //LOG.info("getValueString");
        String str = entry.toString();
        if (str != null) {

            str = str.replace("\n", "<br>");
        }
        return str;
    }
}
