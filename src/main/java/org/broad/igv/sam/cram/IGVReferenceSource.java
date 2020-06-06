/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.sam.cram;


import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import org.apache.log4j.Logger;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.util.ObjectCache;

import java.io.IOException;

/**
 * Provide a reference sequence for CRAM decompression.  Note the rule for MD5 calculation.
 * <p>
 * M5 (sequence MD5 checksum) field of @SQ sequence record in the BAM header is required and UR (URI
 * for the sequence fasta optionally gzipped file) field is strongly advised. The rule for calculating MD5 is
 * to remove any non-base symbols (like \n, sequence name or length and spaces) and upper case the rest.
 */

public class IGVReferenceSource implements CRAMReferenceSource {

    private static Logger log = Logger.getLogger(IGVReferenceSource.class);

    static ObjectCache<String, byte[]> cachedSequences = new ObjectCache<String, byte[]>(2);

    static GenomeChangeListener genomeChangeListener;

    @Override
    public  byte[] getReferenceBases(SAMSequenceRecord record, boolean tryNameVariants) {

        final String name = record.getSequenceName();

        final Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        String chrName = currentGenome.getCanonicalChrName(name);
        Chromosome chromosome = currentGenome.getChromosome(chrName);

        byte[] bases = cachedSequences.get(chrName);

        if (bases == null) {

            try {
                if (bases == null) {
                    if (IGV.hasInstance()) IGV.getInstance().setStatusBarMessage("Loading sequence");
                    bases = currentGenome.getSequence(chrName, 0, chromosome.getLength(), false);

                    // CRAM spec requires upper case
                    for (int i = 0; i < bases.length; i++) {
                        if (bases[i] >= 97) bases[i] -= 32;
                    }
                }

                cachedSequences.put(chrName, bases);
            } finally {
                if (IGV.hasInstance()) IGV.getInstance().setStatusBarMessage("");
            }
        }

        return bases;
    }

    public static class GenomeChangeListener implements IGVEventObserver {

        @Override
        public void receiveEvent(Object event) {
            cachedSequences.clear();
        }
    }

    static {

        genomeChangeListener = new GenomeChangeListener();
        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, genomeChangeListener);
    }


}
