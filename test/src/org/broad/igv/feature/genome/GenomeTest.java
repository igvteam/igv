package org.broad.igv.feature.genome;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class GenomeTest {
    @Test
    public void testGetNCBIName() throws Exception {

        String ncbiID = "gi|125745044|ref|NC_002229.3|";
        String ncbiName = "NC_002229.3";

        assertEquals(ncbiName, Genome.getNCBIName(ncbiID));

    }
}
