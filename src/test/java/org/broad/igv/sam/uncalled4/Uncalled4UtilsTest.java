package org.broad.igv.sam.uncalled4;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMAlignment;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

public class Uncalled4UtilsTest {

    @Test
    public void testComments() throws IOException {

        String bamfile = TestUtils.DATA_DIR + "bam/Uncalled_4/dna_r9_dm_ctl.align.bam";
        String baifile = TestUtils.DATA_DIR + "bam/Uncalled_4/dna_r9_dm_ctl.align.bam.bai";


        String chr = "CM044676.1";
        int start = 0;
        int end = Integer.MAX_VALUE;

        ResourceLocator baiLocator = new ResourceLocator(bamfile);
        baiLocator.setIndexPath(baifile);
        BAMReader baireader = new BAMReader(baiLocator, true);

        CloseableIterator<SAMAlignment> baiiter = baireader.query(chr, start, end, true);

        Map<String, Map<String, Float>> scales = Uncalled4Utils.getScaleFactors(baireader);

        int count = 0;
        while (baiiter.hasNext()) {

            Alignment bamrecord = baiiter.next();


            short[] uc = (short[]) bamrecord.getAttribute("uc");
            short[] ud = (short[]) bamrecord.getAttribute("ud");
            int[] ur = (int[]) bamrecord.getAttribute("ur");

            int nCalls = uc.length;
            assertEquals(nCalls, ud.length);

            int nRanges = ur.length / 2;
            int nBases = 0;
            for (int i = 0; i < nRanges; i++) {
                int s = ur[i * 2];
                int e = ur[i * 2 + 1];
                nBases += e - s;
            }
            assertEquals(nCalls, nBases);


            count++;
        }
        assertTrue("Unexpected data count: " + count, count == 20);

    }

}