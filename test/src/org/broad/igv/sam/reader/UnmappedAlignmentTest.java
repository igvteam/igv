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

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.junit.Test;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jan 25, 2010
 * Time: 10:00:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class UnmappedAlignmentTest {

    @Test
    public void testUnmappedReads() throws Exception {
        File samFile = new File("/Volumes/ifs/seq/picard_aggregation/G2145/TCGA-06-0152-01A-02D/v9/TCGA-06-0152-01A-02D.bam");

        String chr = "chr1";
        int start = 220929208;
        int end = 220929357;

        // Test posA query that includes overlaps (contained == false)
        SAMFileReader reader = new SAMFileReader(samFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        boolean contained = true;
        CloseableIterator<SAMRecord> iter = reader.query(chr, start, end, contained);

        String unmappedName = null;
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            if (record.getReadUnmappedFlag()) {
                System.out.println("Unmapped:     " + record.getReadName());
            }
            if (record.getMateUnmappedFlag() == true) {
                System.out.println("Mate unmapped: " + record.getReadName());
            }
        }
        iter.close();
        reader.close();

        System.out.println("----------------");
        BAMQueryReader reader2 = new BAMQueryReader(samFile);
        CloseableIterator<Alignment> iter2 = reader2.query(chr, start, end, contained);
        unmappedName = null;
        while (iter2.hasNext()) {
            Alignment al = iter2.next();
            if (!al.isMapped()) {
                System.out.println("Unmapped:     " + al.getReadName());
            }
            if (!al.getMate().isMapped()) {
                System.out.println("Mate unmapped: " + al.getReadName());
            }
        }
        iter2.close();
        reader2.close();

    }
}
