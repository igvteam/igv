package org.igv.alignment.fiberseq;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.igv.alignment.SAMAlignment;
import org.igv.alignment.fiberseq.FiberseqAnnotations.Interval;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * Fiber-seq annotations read from the fibertools test BAMs in test/data/bam/fiberseq.  Expected intervals were
 * computed independently of IGV by walking each read's CIGAR.  Intervals are listed in molecular order, which is
 * descending reference order on reverse-strand reads.
 */
public class FiberseqBamTest {

    private static final String DIR = TestUtils.DATA_DIR + "bam/fiberseq/";

    // The same reverse-strand read is in msp_nuc.bam (legacy tags) and ma_spelled.bam (Ma tags)
    private static final String MSP_NUC_READ = "m54329U_220210_004342/140313102/ccs";

    /**
     * Annotations for every mapped record with fiber-seq tags, keyed by read name.
     */
    private static Map<String, FiberseqAnnotations> load(String file) throws IOException {
        Map<String, FiberseqAnnotations> annotations = new LinkedHashMap<>();
        try (SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(new File(DIR + file))) {
            for (SAMRecord record : reader) {
                if (!record.getReadUnmappedFlag() && FiberseqAnnotations.hasTags(record::getAttribute)) {
                    annotations.put(record.getReadName(), new SAMAlignment(record).getFiberseqAnnotations());
                }
            }
        }
        return annotations;
    }

    private static List<Interval> fire(FiberseqAnnotations a) {
        return a.getMsps().stream().filter(i -> i.quality() > 0).toList();
    }

    @Test
    public void legacyTags() throws IOException {
        FiberseqAnnotations a = load("msp_nuc.bam").get(MSP_NUC_READ);
        assertEquals(35, a.getNucleosomes().size());
        assertEquals(36, a.getMsps().size());
        assertEquals(List.of(new Interval(3081870, 3082038, 0), new Interval(3081715, 3081853, 0)),
                a.getNucleosomes().subList(0, 2));
        assertEquals(List.of(new Interval(3082038, 3082108, 0), new Interval(3081853, 3081870, 0)),
                a.getMsps().subList(0, 2));
    }

    @Test
    public void maTagsMatchLegacyTagsForTheSameRead() throws IOException {
        FiberseqAnnotations legacy = load("msp_nuc.bam").get(MSP_NUC_READ);
        FiberseqAnnotations ma = load("ma_spelled.bam").get(MSP_NUC_READ);
        assertEquals(legacy.getNucleosomes(), ma.getNucleosomes());
        assertEquals(legacy.getMsps(), ma.getMsps());
    }

    @Test
    public void legacyFireQualities() throws IOException {
        FiberseqAnnotations a = load("nuc_example.bam").get("d203cb1f-3b46-4670-8576-19ce6f74cdae");
        assertEquals(18, a.getNucleosomes().size());
        assertEquals(17, a.getMsps().size());
        assertEquals(new Interval(3370742, 3370852, 0), a.getNucleosomes().get(0));
        assertEquals(List.of(new Interval(3370494, 3370742, 245), new Interval(3369724, 3370018, 245)), fire(a));
    }

    @Test
    public void uppercaseMaTagsWithFireQualities() throws IOException {
        Map<String, FiberseqAnnotations> reads = load("NAPA.bam");
        assertEquals(154, reads.size());
        int nucleosomes = 0, msps = 0, fires = 0;
        for (FiberseqAnnotations a : reads.values()) {
            if (a != null) {
                nucleosomes += a.getNucleosomes().size();
                msps += a.getMsps().size();
                fires += fire(a).size();
            }
        }
        assertEquals(12721, nucleosomes);
        assertEquals(12765, msps);
        assertEquals(397, fires);

        FiberseqAnnotations forward = reads.get("m54329U_210814_130637/54723395/ccs");
        assertEquals(152, forward.getNucleosomes().size());
        assertEquals(152, forward.getMsps().size());
        assertEquals(new Interval(47480236, 47480446, 0), forward.getNucleosomes().get(0));
        assertEquals(new Interval(47480235, 47480236, 0), forward.getMsps().get(0));
        assertEquals(9, fire(forward).size());
        assertEquals(new Interval(47481337, 47481485, 243), fire(forward).get(0));

        FiberseqAnnotations reverse = reads.get("m84039_230404_003541_s3/70845505/ccs");
        assertEquals(170, reverse.getNucleosomes().size());
        assertEquals(170, reverse.getMsps().size());
        assertEquals(new Interval(47518555, 47518659, 0), reverse.getNucleosomes().get(0));
        assertEquals(new Interval(47518537, 47518555, 0), reverse.getMsps().get(0));
        assertEquals(2, fire(reverse).size());
        assertEquals(new Interval(47518265, 47518421, 250), fire(reverse).get(0));
    }
}
