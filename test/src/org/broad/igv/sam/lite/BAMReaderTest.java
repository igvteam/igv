package org.broad.igv.sam.lite;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.*;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Created by jrobinso on 3/14/17.
 */
public class BAMReaderTest {

    @Test
    public void readAlignments() throws Exception {

        String bamPath = TestUtils.DATA_DIR + "bam/four.reads.bam";
        String chr = "13";
        int beg = 32353128;
        int end = 32353284;

        BAMReader bamReader = new BAMReader(bamPath);

        List<Alignment> alignmentContainer = bamReader.readAlignments(chr, beg, end);

        String attr = ((BAMAlignment) alignmentContainer.get(0)).getAttributeString(false);

        assertEquals(4, alignmentContainer.size());

    }

    @Test
    public void iterateAlignments() throws Exception {

        String bamPath = TestUtils.DATA_DIR + "bam/four.reads.bam";
        String chr = "13";
        int beg = 32353128;
        int end = 32353284;

        BAMReader bamReader = new BAMReader(bamPath);

        CloseableIterator<Alignment> iter = bamReader.query(chr, beg, end, false);

        List<Alignment> alignmentContainer = new ArrayList<>();
        while(iter.hasNext()) {
            alignmentContainer.add(iter.next());
        }

        String attr = ((BAMAlignment) alignmentContainer.get(0)).getAttributeString(false);

        assertEquals(4, alignmentContainer.size());

    }

    @Test
    public void readHeader() throws Exception {

        String bamPath = TestUtils.DATA_DIR + "bam/gstt1_sample.bam";

        BAMReader bamReader = new BAMReader(bamPath);

        bamReader.readHeader();

        assertEquals(0, (int) bamReader.chrToIndex.get("chr1"));


    }
}
