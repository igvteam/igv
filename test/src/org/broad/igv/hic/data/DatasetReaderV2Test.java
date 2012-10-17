package org.broad.igv.hic.data;

import org.junit.Test;

/**
 * @author jrobinso
 *         Date: 10/17/12
 *         Time: 10:36 AM
 */
public class DatasetReaderV2Test {
    @Test
    public void testRead() throws Exception {
        String path = "/Users/jrobinso/projects/hic/New/MiSeq_GM_HINDIII.hic";

        DatasetReaderV2 reader = new DatasetReaderV2(path);

        Dataset ds = reader.read();

    }
}
