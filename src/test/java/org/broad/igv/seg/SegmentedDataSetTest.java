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

public class SegmentedDataSetTest {

    static final String testFile = "test/data/seg/Broad.080528.subtypes.seg.gz";
    private static SegmentedDataSet segmentedDataSet;

    @BeforeClass
    public static void setup() throws IOException {
        Genome genome = TestUtils.mockUCSCGenome();
        ResourceLocator locator = new ResourceLocator(testFile);
        segmentedDataSet = SegmentFileParser.loadSegments(locator, genome);
    }


    @Test
    public void getChromosomes() {
        Set<String> chromosomes = segmentedDataSet.getChromosomes();
        Assert.assertTrue(chromosomes.size() > 0);
    }

    @Test
    public void getSegments() {

    }

    @Test
    public void getSampleNames() {
        final List<String> sampleNames = segmentedDataSet.getSampleNames();
        Assert.assertTrue(sampleNames.size() > 0);
    }

    @Test
    public void isLogNormalized() {
        Assert.assertTrue(segmentedDataSet.isLogNormalized());
    }

    @Test
    public void getDataMax() {
    }

    @Test
    public void getDataMin() {
    }

    @Test
    public void getWholeGenomeScores() {
        final List<String> sampleNames = segmentedDataSet.getSampleNames();
        for (String sampleName : sampleNames) {
            final List<LocusScore> wholeGenomeScores = segmentedDataSet.getWholeGenomeScores(sampleName);
            Assert.assertTrue(wholeGenomeScores.size() > 0);
        }
    }
}