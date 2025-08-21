package org.broad.igv.htsget;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import junit.framework.Assert;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMAlignment;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class HtsgetBAMReaderTest {


    /**
     * Minimal test, just verifies that some bam records are returned.
     *
     * @throws Exception
     */
    @Test
    public void testQueryAlignments() throws Exception {

        String url = "https://htsget.ga4gh-demo.org/reads/htsnexus_test_NA12878";
        String chr = "11";
        int start = 5020134;
        int end = 5020614;

        ResourceLocator locator = new ResourceLocator(url);
        locator.setHtsget(true);

        BAMReader bamreader = new BAMReader(locator, false);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, true);

        int count = 0;
        List<Alignment> alignmentList = new ArrayList<>();
        while (bamiter.hasNext()) {
            Alignment bamrecord = bamiter.next();
            if (bamrecord.getEnd() > start && bamrecord.getStart() < end &&
                    bamrecord.isVendorFailedRead() == false &&
                   bamrecord.isMapped()
            ) {
                alignmentList.add(bamrecord);
            }
            count++;
        }

        assertTrue("No data retrieved", alignmentList.size() > 0);
    }


//    public static void main(String [] args) throws Exception {
//        (new HtsgetBAMReaderTest()).testQueryAlignments();
//    }


}