package org.broad.igv.track;

import htsjdk.samtools.seekablestream.SeekableStream;
import junit.framework.TestCase;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

public class FileFormatUtilsTest extends TestCase {

    public void testIsBAM() throws Exception {
        String path = "https://1000genomes.s3.amazonaws.com/phase3/data/HG01879/exome_alignment/HG01879.mapped.ILLUMINA.bwa.ACB.exome.20120522.bam";
        boolean b = FileFormatUtils.isBAM(path);
        assertTrue(b);

    }

    public void testDetermineFormat() throws Exception {
        String bamFIle = TestUtils.DATA_DIR + "bam/NA12878.SLX.sample.bam";
        String format = FileFormatUtils.determineFormat(bamFIle);
        assertEquals("bam", format);

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_bai_index.cram";
        format = FileFormatUtils.determineFormat(cramFile);
        assertEquals("cram", format);

        String vcfFile = TestUtils.DATA_DIR + "vcf/ex2.vcf";
        format = FileFormatUtils.determineFormat(vcfFile);
        assertEquals("vcf", format);

        String gffFile = TestUtils.DATA_DIR + "gff/gene.sorted.gff3";
        format = FileFormatUtils.determineFormat(gffFile);
        assertEquals("gff3", format);

        String tdfFile = TestUtils.DATA_DIR + "tdf/NA12878.SLX.egfr.sam.tdf";
        format = FileFormatUtils.determineFormat(tdfFile);
        assertEquals("tdf", format);

        String unknown = TestUtils.DATA_DIR + "cram/cram_with_crai_index.cram.crai";
        format = FileFormatUtils.determineFormat(unknown);
        assertNull(format);
    }
}