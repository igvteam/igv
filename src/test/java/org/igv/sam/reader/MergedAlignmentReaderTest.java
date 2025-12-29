package org.igv.sam.reader;

import htsjdk.samtools.util.CloseableIterator;
import org.igv.AbstractHeadlessTest;
import org.igv.sam.Alignment;
import org.igv.tools.IGVToolsTest;
import org.igv.util.FileUtils;
import org.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jacob
 * @since 2012/01/25
 */
public class MergedAlignmentReaderTest extends AbstractHeadlessTest {

    public static String[] generateRepLargebamsList(File listFile) throws IOException {
        listFile.delete();
        listFile.deleteOnExit();
        String listPath = listFile.getPath();

        return IGVToolsTest.generateRepLargebamsList(listPath, "HG00171.hg18.bam", 2);
    }

    @Test
    @Ignore("Requires largedata bundle")
    public void testSimpleRead() throws Exception {

        File listFile = new File(TestUtils.LARGE_DATA_DIR, "2largebams.bam.list");

        String[] actfiles = generateRepLargebamsList(listFile);
        int start = 151667156;
        int end = start + 10000;
        int num_combined = 0;
        int num_sep = 0;
        Alignment align;

        AlignmentReader mergedReader = AlignmentReaderFactory.getBamListReader(listFile.getAbsolutePath(), false);
        CloseableIterator<Alignment> combData = mergedReader.query("chr1", start, end, false);

        Map<Float, Integer> combinedCounts = new HashMap();
        while (combData.hasNext()) {
            align = combData.next();
            num_combined += 1;
            float score = align.getScore();
            Integer oldcount = combinedCounts.get(score);
            int newcount = oldcount == null ? 1 : oldcount + 1;
            combinedCounts.put(score, newcount);
        }
        mergedReader.close();

        BufferedReader in = new BufferedReader(new FileReader(listFile));
        String singfile = in.readLine();
        singfile = FileUtils.getAbsolutePath(singfile, listFile.getPath());
        AlignmentReader singReader = AlignmentReaderFactory.getReader(singfile, false);

        CloseableIterator<Alignment> singData = singReader.query("chr1", start, end, false);

        Map<Float, Integer> singCounts = new HashMap();
        while (singData.hasNext()) {
            align = singData.next();
            num_sep += 1;
            float score = align.getScore();
            Integer oldcount = singCounts.get(score);
            int newcount = oldcount == null ? 1 : oldcount + 1;
            singCounts.put(score, newcount);
        }

        singReader.close();

        assertEquals(2 * num_sep, num_combined);
        assertEquals(singCounts.size(), combinedCounts.size());
        for (Float score : singCounts.keySet()) {
            assertEquals(2 * singCounts.get(score), (int) combinedCounts.get(score));
        }
    }


    @Test
    public void testSortOrder() throws Exception {
        String listPath = TestUtils.DATA_DIR + "bam/test.unindexed.bam.list";

        AlignmentReader mergedReader = AlignmentReaderFactory.getBamListReader(listPath, false);
        CloseableIterator<Alignment> iter = mergedReader.iterator();

        int lastPosition = 0;
        String lastChr = "";
        while (iter.hasNext()) {
            Alignment a = iter.next();
            String chr = a.getChr();
            int pos = a.getAlignmentStart();
            assertTrue(chr.compareTo(lastChr) >= 0);
            if(lastChr.equals(chr)) {
                assertTrue(pos >= lastPosition);
            }
            lastChr = chr;
            lastPosition = pos;
        }
        iter.close();
    }

}
