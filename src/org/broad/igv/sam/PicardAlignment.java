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

package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.color.ColorUtilities;

import java.util.List;

/**
 * @author jrobinso
 */
public class PicardAlignment extends SAMAlignment implements Alignment {

    private static Logger log = Logger.getLogger(PicardAlignment.class);

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


    private int flags;

    /**
     * Picard object upon which this PicardAlignment is based
     */
    private SAMRecord record;

    public PicardAlignment(SAMRecord record) {
        super();

        this.record = record;

        this.flags = record.getFlags();

        String refName = record.getReferenceName();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.chr = genome == null ? refName : genome.getChromosomeAlias(refName);

        // SAMRecord is 1 based inclusive.  IGV is 0 based exclusive.

        this.end = record.getAlignmentEnd();   // might be modified later for soft clipping
        this.start = record.getAlignmentStart() - 1;   // might be modified later for soft clipping


        if (record.getReadPairedFlag()) {
            String mateReferenceName = record.getMateReferenceName();
            String mateChr = genome == null ? mateReferenceName : genome.getChromosomeAlias(mateReferenceName);
            this.setMate(new ReadMate(mateChr,
                    record.getMateAlignmentStart() - 1,
                    record.getMateNegativeStrandFlag(),
                    record.getMateUnmappedFlag()));
        }

        String keySequence = null;
        SAMFileHeader header = record.getHeader();
        String flowOrder = null;
        if (header != null) {
            readGroup = (String) record.getAttribute("RG");
            if (readGroup != null) {
                SAMReadGroupRecord rgRec = header.getReadGroup(readGroup);
                if (rgRec != null) {
                    this.sample = rgRec.getSample();
                    this.library = rgRec.getLibrary();
                    flowOrder = rgRec.getFlowOrder();
                    keySequence = rgRec.getKeySequence();
                }
            }
        }

        Object colorTag = record.getAttribute("YC");
        if (colorTag != null) {
            try {
                color = ColorUtilities.stringToColor(colorTag.toString(), null);
            } catch (Exception e) {
                log.error("Error interpreting color tag: " + colorTag, e);
            }
        }

        setPairOrientation();
        setPairStrands();
        createAlignmentBlocks(record.getCigarString(), record.getReadBases(), record.getBaseQualities(),
                getFlowSignals(flowOrder, keySequence), flowOrder, this.getFlowSignalsStart());


    }      // End constructor

    /**
     * @return The SAMRecord which created this PicardAlignment
     */
    public SAMRecord getRecord() {
        return this.record;
    }

    public Object getAttribute(String key) {
        // SAM alignment tag keys must be of length 2
        return key.length() == 2 ? record.getAttribute(key) :
                (key.equals("TEMPLATE_ORIENTATION") ? pairOrientation : null);
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

    @Override
    public boolean isSupplementary() {
        return (flags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }


    public boolean isVendorFailedRead() {
        return (flags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }
    @Override
    public boolean isPrimary() {
        return (flags & NOT_PRIMARY_ALIGNMENT_FLAG) == 0;
    }

    @Override
    public String toString() {
        return record.getSAMString();
    }

    @Override
    public String getReadName() {
        return record.getReadName();
    }

    @Override
    public int getMappingQuality() {
        return record.getMappingQuality();
    }

    @Override
    public int getInferredInsertSize() {
        return record.getInferredInsertSize();
    }

    @Override
    public String getCigarString() {
        return record.getCigarString();
    }

    @Override
    public int getReadLength() {
        return record.getReadString().length();
    }

    @Override
    public String getReadSequence() {
        return record.getReadString();
    }

    @Override
    public int getAlignmentStart() {
        return record.getAlignmentStart() - 1;
    }

    @Override
    public int getAlignmentEnd() {
        return record.getAlignmentEnd();
    }

    protected String getAttributeString(boolean truncate) {

        StringBuffer buf = new StringBuffer();
        SAMRecord record = getRecord();
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
        return buf.toString();
    }

    /**
     * @param flowOrder   the flow order corresponding to this read
     * @param keySequence sequence the key sequence corresponding to this read
     * @return the flow signals in 100x format (SFF), only if they exist (FZ tag),
     * if the key sequence and flow order are found in the read group header tag
     * (RG.KS and RG.FO).  Note: the array proceeds in the sequencing direction.
     */
    public short[] getFlowSignals(String flowOrder, String keySequence) {
        short[] r = null;
        int i;
        int startFlow, keySignalOverlap;
        char firstBase;

        if (null == flowOrder || null == keySequence) {
            return null;
        }

        startFlow = this.getFlowSignalsStart();
        if (startFlow < 0) {
            return null;
        }

        // get the # of bases that the first base in the read overlaps with the last base(s) in the key
        SAMRecord record = getRecord();
        if (this.isNegativeStrand()) {
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
