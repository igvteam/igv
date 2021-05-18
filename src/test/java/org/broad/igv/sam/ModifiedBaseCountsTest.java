package org.broad.igv.sam;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.*;

public class ModifiedBaseCountsTest {

    @Test
    public void incrementCounts() throws IOException {

        String bamfile = TestUtils.DATA_DIR + "bam/chr20_mod_call_sample.bam";
        String chr = "20";
        int start = 13846149;
        int end = 13846294;

        BAMReader bamreader = new BAMReader(new ResourceLocator(bamfile), true);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, false);
        int readCount = 0;

        ModifiedBaseCounts counts = new ModifiedBaseCounts();
        while (bamiter.hasNext()) {
            Alignment alignment = bamiter.next();
            counts.incrementCounts(alignment);
            readCount++;
        }
        assertTrue("No data retrieved:  " + readCount, readCount > 0);

        int[] expectedPositions = {13846181, 13846182, 13846227, 13846228, 13846232, 13846233, 13846234};
        int[] expectedCounts =    {1,        1,        1,        2,        1,        1,        1       };

        for(int i=0; i<expectedPositions.length; i++) {
            int c = counts.getCount(expectedPositions[i] - 1, "m");
            assertEquals("Unexpected count at position " + expectedPositions[i], expectedCounts[i], c);
        }

       // counts.dump();

    }
}