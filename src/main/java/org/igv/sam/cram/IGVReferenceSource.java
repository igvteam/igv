package org.igv.sam.cram;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import org.igv.logging.*;
import org.igv.event.GenomeChangeEvent;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.feature.Chromosome;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.ui.IGV;
import org.igv.util.ObjectCache;

import java.util.Arrays;
import java.util.HashMap;

/**
 * Provide a reference sequence for CRAM decompression.   The rule for calculating MD5 is
 * to remove any non-base symbols (like \n, sequence name or length and spaces) and upper case the rest.
 */

public class IGVReferenceSource implements CRAMReferenceSource {

    private static Logger log = LogManager.getLogger(IGVReferenceSource.class);

    @Override
    public byte[] getReferenceBases(SAMSequenceRecord record, boolean tryNameVariants) {
        return getReferenceBasesByRegion(record, 0, record.getSequenceLength());
    }

    @Override
    public byte[] getReferenceBasesByRegion(final SAMSequenceRecord sequenceRecord, final int zeroBasedStart, final int requestedRegionLength) {
        final String name = sequenceRecord.getSequenceName();
        final Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        String chrName = currentGenome.getCanonicalChrName(name);
        byte[] bases = currentGenome.getSequence(chrName, zeroBasedStart, zeroBasedStart + requestedRegionLength, true);

        if (bases != null) {
            // CRAM spec requires upper case
            for (int i = 0; i < bases.length; i++) {
                if (bases[i] >= 97) bases[i] -= 32;
            }
        }
        return bases;
    }

}
