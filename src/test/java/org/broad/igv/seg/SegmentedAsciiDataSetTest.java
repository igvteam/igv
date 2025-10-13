package org.broad.igv.seg;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.Set;

public class SegmentedAsciiDataSetTest {

    static final String testFile = "test/data/seg/Broad.080528.subtypes.seg.gz";
    private static SegmentedAsciiDataSet segmentedAsciiDataSet;

    @BeforeClass
    public static void setup() throws IOException {
        Genome genome = TestUtils.mockUCSCGenome();
        ResourceLocator locator = new ResourceLocator(testFile);
        segmentedAsciiDataSet = SegmentFileParser.loadSegments(locator, genome);
    }


    @Test
    public void getChromosomes() {
        Set<String> chromosomes = segmentedAsciiDataSet.getChromosomes();
        Assert.assertTrue(chromosomes.size() > 0);
    }

    @Test
    public void getSegments() {

    }

    @Test
    public void getSampleNames() {
        final List<String> sampleNames = segmentedAsciiDataSet.getSampleNames();
        Assert.assertTrue(sampleNames.size() > 0);
    }

    @Test
    public void isLogNormalized() {
        Assert.assertTrue(segmentedAsciiDataSet.isLogNormalized());
    }

    @Test
    public void getDataMax() {
    }

    @Test
    public void getDataMin() {
    }

    @Test
    public void getWholeGenomeScores() {
        final List<String> sampleNames = segmentedAsciiDataSet.getSampleNames();
        for (String sampleName : sampleNames) {
            final List<LocusScore> wholeGenomeScores = segmentedAsciiDataSet.getWholeGenomeScores(sampleName);
            Assert.assertTrue(wholeGenomeScores.size() > 0);
        }
    }
}