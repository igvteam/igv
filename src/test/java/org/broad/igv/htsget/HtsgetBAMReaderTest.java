package org.broad.igv.htsget;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import junit.framework.Assert;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMAlignment;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.Iterator;

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

        String endpoint = "htsget://htsget.ga4gh.org/reads/giab.NA12878.NIST7086.1";
        String chr = "chr8";
        int start = 128734098;
        int end = 128763217;

        ResourceLocator locator = new ResourceLocator(endpoint);
        locator.setHtsget(true);

        BAMReader bamreader = new BAMReader(locator, false);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, true);

        int count = 0;
        while (bamiter.hasNext()) {
            Alignment bamrecord = bamiter.next();
            count++;
        }
        assertTrue("No data retrieved", count > 0);
    }


//    public static void main(String [] args) throws Exception {
//        (new HtsgetBAMReaderTest()).testQueryAlignments();
//    }


}