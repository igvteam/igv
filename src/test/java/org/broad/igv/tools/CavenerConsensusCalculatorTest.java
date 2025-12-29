package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.DenseAlignmentCounts;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;

import static junit.framework.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Jun-24
 */
@RunWith(Parameterized.class)
public class CavenerConsensusCalculatorTest extends AbstractHeadlessTest {

    private String inputBases;
    private char expectedConsensus;

    public CavenerConsensusCalculatorTest(String inputBases, char expectedConsensus) {
        this.inputBases = inputBases;
        this.expectedConsensus = expectedConsensus;

        //System.out.println("Input: " + this.inputBases);
    }

    @Test
    public void testCalcCorrectly() throws Exception {
        CavenerConsensusCalculator calc = new CavenerConsensusCalculator();
        AlignmentCounts counts = createCounts();
        assertEquals(this.expectedConsensus, calc.calculateConsensusBase(counts, 0));
    }

    private AlignmentCounts createCounts() {
        MockAlignmentCounts counts = new MockAlignmentCounts(0, 1);
        for (char c : this.inputBases.toCharArray()) {
            counts.incPositionCount(0, (byte) c);
        }
        counts.finish();

        return counts;
    }


    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        Object[][] data = new Object[][]{
                {"aaa", 'a'}, {"ccccccg", 'c'}, {"g", 'g'}, {"n", 'n'},
                {"aac", 'a'}, {"acgt", 'n'},
                {"aacc", 'm'}, {"aacct", 'm'}, {"ggaaactgga", 'r'}, {"ggggccccatat", 'n'}


        };
        return Arrays.asList(data);
    }

    private static class MockAlignmentCounts extends DenseAlignmentCounts {

        public MockAlignmentCounts(int start, int end) {
            super(start, end, null);
        }

        public void incPositionCount(int pos, byte b) {
            super.incPositionCount(pos, b, (byte) 0, false);
        }
    }
}
