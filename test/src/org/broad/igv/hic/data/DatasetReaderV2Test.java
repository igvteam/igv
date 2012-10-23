/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.hic.data;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.ChromosomeImpl;
import org.junit.Ignore;
import org.junit.Test;

/**
 * @author jrobinso
 *         Date: 10/17/12
 *         Time: 10:36 AM
 */
@Ignore
public class DatasetReaderV2Test {

    @Test
    public void testRead() throws Exception {
        String path = "/Users/jrobinso/projects/hic/New/MiSeq_GM_HINDIII.hic";

        DatasetReaderV2 reader = new DatasetReaderV2(path);

        Dataset ds = reader.read();

        Chromosome chr1 = new ChromosomeImpl(14, "14", 0);
        Chromosome chr2 = chr1;

        Matrix matrix = ds.getMatrix(chr1, chr2);

        MatrixZoomData mzd = matrix.getObservedMatrix(0);

        Block b = mzd.getBlock(0);

        b.getContactRecords();

    }
}
