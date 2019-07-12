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

package org.broad.igv.sam.cram;

import htsjdk.samtools.*;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.stream.IGVSeekableBufferedStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Created by jrobinso on 5/25/16.
 */

public class IGVReferenceSourceTest {

    @Before
    public void setUp() throws Exception {

    }

    @Test
    public void testGetReferenceBases() throws Exception {

        String fastaURL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa";

        GenomeManager.getInstance().loadGenome(fastaURL, null);

        IGVReferenceSource refSource = new IGVReferenceSource();
        SAMSequenceRecord rec = new SAMSequenceRecord("22", 50818468);
        byte[] bases = refSource.getReferenceBases(rec, false);

        assertEquals(50818468, bases.length);

        assertEquals('N', bases[0]);
        assertEquals('G', bases[27198882]);
    }


    @Test
    public void testGetReferenceBasesCompressed() throws Exception {

        String fastaURL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.gz";

        GenomeManager.getInstance().loadGenome(fastaURL, null);

        IGVReferenceSource refSource = new IGVReferenceSource();
        SAMSequenceRecord rec = new SAMSequenceRecord("22", 50818468);
        byte[] bases = refSource.getReferenceBases(rec, false);

        assertEquals(50818468, bases.length);

        assertEquals('N', bases[0]);
        assertEquals('G', bases[27198882]);
    }

//    @Test
//    public void testCompressedTiming() throws Exception {
//
//        String fastaURL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa";
//        GenomeManager.getInstance().loadGenome(fastaURL, null);
//        IGVReferenceSource refSource = new IGVReferenceSource();
//        SAMSequenceRecord rec = new SAMSequenceRecord("1", 248956422);
//
//        long t = System.currentTimeMillis();
//        byte[] bases = refSource.getReferenceBases(rec, false);
//        assertEquals(248956422, bases.length);
//        long dt = System.currentTimeMillis() - t;
//
//
//        fastaURL = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.gz";
//        GenomeManager.getInstance().loadGenome(fastaURL, null);
//        refSource = new IGVReferenceSource();
//        rec = new SAMSequenceRecord("1", 248956422);
//
//        long t1 = System.currentTimeMillis();
//        bases = refSource.getReferenceBases(rec, false);
//        assertEquals(248956422, bases.length);
//        long dt1 = System.currentTimeMillis() - t1;
//
//        System.out.println(dt + "    " + dt1);
//
//        assertTrue(dt1 < dt);
//
//    }
}
