package org.broad.igv.sam.mods;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMAlignment;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

public class BaseModificationCountsTest {

    @Before
    public void setup() {
        PreferencesManager.getPreferences().put(Constants.BASEMOD_VALIDATE_BASE_COUNT, "false");
    }

    @Test
    public void incrementCounts() throws IOException {

        String bamfile = TestUtils.DATA_DIR + "bam/chr20_mod_call_sample.bam";
        String chr = "20";
        int start = 13846149;
        int end = 13846294;

        BAMReader bamreader = new BAMReader(new ResourceLocator(bamfile), true);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, false);
        int readCount = 0;

        BaseModificationCounts counts = new BaseModificationCounts();
        while (bamiter.hasNext()) {
            Alignment alignment = bamiter.next();
            counts.incrementCounts(alignment);
            readCount++;
        }
        assertTrue("No data retrieved:  " + readCount, readCount > 0);

        int[] expectedPositions = {13846181, 13846182, 13846227, 13846228, 13846232, 13846233, 13846234};
        int[] expectedCounts =    {1,        1,        1,        2,        1,        1,        1       };

        BaseModificationKey key =  BaseModificationKey.getKey('C', '+', "m");

        boolean includeNoMods = false;
        for(int i=0; i<expectedPositions.length; i++) {
            int c = counts.getCount(expectedPositions[i] - 1, key, 0, includeNoMods);
            assertEquals("Unexpected count at position " + expectedPositions[i], expectedCounts[i], c);
        }
    }

    @Test
    public void incrementCounts2() throws IOException {

        String bamfile = "https://www.dropbox.com/s/q32hk7tvsejryjt/HG002_chr11_119076212_119102218_2.bam";
        String indexFile = "https://www.dropbox.com/s/ax1ljny7ja5fdcu/HG002_chr11_119076212_119102218_2.bam.bai";
        String chr = "chr11";
        int start = 119094722;
        int end = 119094724;
        boolean includeNoMods = true;

        ResourceLocator locator = new ResourceLocator(bamfile);
        locator.setIndexPath(indexFile);
        BAMReader bamreader = new BAMReader(locator, true);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, false);
        int readCount = 0;

        BaseModificationCounts counts = new BaseModificationCounts();
        while (bamiter.hasNext()) {
            Alignment alignment = bamiter.next();
            counts.incrementCounts(alignment);
            readCount++;
        }
        assertTrue("No data retrieved:  " + readCount, readCount > 0);

        BaseModificationKey cmKey = BaseModificationKey.getKey('C', '+', "m");
        int aboveThreshold = counts.getCount(119094723, cmKey, 0.5f, includeNoMods);
        assertEquals("Counts above threshold", 3, aboveThreshold);

        BaseModificationKey noModKey = BaseModificationKey.getKey('C', '+', "NONE_C");
        int belowThreshold = counts.getCount( 119094723, noModKey, 0.5f, includeNoMods);
        assertEquals("Counts below threshold", 12, belowThreshold);
    }
}