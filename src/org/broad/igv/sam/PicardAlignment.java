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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.color.ColorUtilities;

import java.util.List;

/**
 * @author jrobinso
 */
public class PicardAlignment extends SAMAlignment implements Alignment {

    private static Logger log = Logger.getLogger(PicardAlignment.class);

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

        this.alignmentStart = record.getAlignmentStart() - 1;
        this.alignmentEnd = Math.max(alignmentStart, record.getAlignmentEnd());
        this.end = alignmentEnd;   // might be modified later for soft clipping
        this.start = this.alignmentStart;   // might be modified later for soft clipping
        this.mappingQuality = record.getMappingQuality();
        this.readName = record.getReadName().trim();
        this.inferredInsertSize = record.getInferredInsertSize();

        this.cigarString = record.getCigarString();
        this.readSequence = record.getReadString();

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
                color = ColorUtilities.stringToColor(colorTag.toString());
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


    @Override
    public String toString() {
        return record.getSAMString();
    }

    public String getClipboardString(double location) {
        return getValueStringImpl(location, false);
    }

    String getValueStringImpl(double position, boolean truncate) {

        StringBuffer buf = new StringBuffer(super.getValueString(position, null));
        SAMRecord record = getRecord();
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
                buf.append("<br>FAILED Vendor quality check");
                sectionBreak = true;
            }
            if (sectionBreak) {
                buf.append("<br>-------------------");
            }
        }

        if(record.getSupplementaryAlignmentFlag()){
            buf.append("<br>Supplementary alignment (chimeric)");
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
