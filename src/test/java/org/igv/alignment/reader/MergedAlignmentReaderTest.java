package org.igv.alignment.reader;

import htsjdk.samtools.util.CloseableIterator;
import org.igv.AbstractHeadlessTest;
import org.igv.alignment.Alignment;
import org.igv.util.TestUtils;
import org.junit.Test;


import static junit.framework.Assert.assertTrue;

/**
 * @author jacob
 * @since 2012/01/25
 */
public class MergedAlignmentReaderTest extends AbstractHeadlessTest {


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
