package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.junit.Test;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 2/13/13
 *         Time: 2:21 PM
 */
public class ChromosomeComparatorTest {


    // Test chr names which resolve to numbers too large for "int".  In response to bug report.

    @Test
    public void testCompareLargeNumbers() throws Exception {
        long num1 = 7000000037415152l;
        long num2 = 7000000037415153l;
        ChromosomeComparator comp = new ChromosomeComparator(0);
        Chromosome chr1 = new Chromosome(0, String.valueOf(num1), 1 );
        Chromosome chr2 = new Chromosome(0, String.valueOf(num2), 1 );
        int value = comp.compare(chr1, chr2);
        assertTrue(value < 0);
    }
}
