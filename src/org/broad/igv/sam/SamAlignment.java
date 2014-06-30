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

/**
 * @author jrobinso
 */
public class SamAlignment extends AbstractAlignment implements Alignment {

    private static Logger log = Logger.getLogger(SamAlignment.class);

    public SamAlignment(SAMRecord record) {
        super();
        String keySequence = null;

        this.negativeStrand = record.getReadNegativeStrandFlag();
        this.duplicate = record.getDuplicateReadFlag();
        this.mapped =  !record.getReadUnmappedFlag();
        this.readLength = record.getReadLength();
        this.paired = record.getReadPairedFlag();
        this.properPair = record.getProperPairFlag();
        this.firstOfPair = record.getFirstOfPairFlag();
        this.secondOfPair = record.getSecondOfPairFlag();
        this.cigarString = record.getCigarString();
        this.readSequence = record.getReadString();
        this.primary =   !record.getNotPrimaryAlignmentFlag();
        this.supplementary = record.getSupplementaryAlignmentFlag();

        String refName = record.getReferenceName();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.chr = genome == null ? refName : genome.getChromosomeAlias(refName);

        // SAMRecord is 1 based inclusive.  IGV is 0 based exclusive.
        this.start = this.alignmentStart;   // might be modified later for soft clipping
        this.end = alignmentEnd;   // might be modified later for soft clipping
        this.setMappingQuality(record.getMappingQuality());
        this.readName = record.getReadName().trim();
        this.setInferredInsertSize(record.getInferredInsertSize());

        setMatePair(genome);
        setPairOrientation();
        setPairStrands();

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

        createAlignmentBlocks(record.getCigarString(), record.getReadBases(), record.getBaseQualities(), decodeReduceCounts(record),
                getFlowSignals(flowOrder, keySequence), flowOrder, this.getFlowSignalsStart());

        Object colorTag = record.getAttribute("YC");
        if (colorTag != null) {
            try {
                color = ColorUtilities.stringToColor(colorTag.toString());
            } catch (Exception e) {
                log.error("Error interpreting color tag: " + colorTag, e);
            }
        }
    }      // End constructor


}
