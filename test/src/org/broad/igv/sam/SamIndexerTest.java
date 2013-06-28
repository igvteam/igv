/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.AlignmentIndexer;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.sam.reader.SamIndexer;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class SamIndexerTest extends AbstractHeadlessTest{

    /**
     * Test of createSamIndex method, of class SamIndexer.
     */
    @Test
    public void createSamIndex() throws Exception {
        String samFilePath = TestUtils.DATA_DIR + "sam/test_2.sam";
        File samFile = new File(samFilePath);
        File idxFile = new File(samFilePath + ".sai");
        idxFile.deleteOnExit();
        int tileWidth = 100;
        AlignmentIndexer instance = new SamIndexer(samFile, null, null);
        //AlignmentIndexer.getInstance(samFile, null, null);

        FeatureIndex result = instance.createSamIndex(idxFile, tileWidth);
        result.store(idxFile);

        SAMReader reader = new SAMReader(samFilePath, true);

        int start = 125963669;
        CloseableIterator<Alignment> iter = reader.query("chr3", start, start + 200, false);
        int count = 0;

        while (iter.hasNext()) {
            Alignment rec = iter.next();
            count++;
        }

        iter.close();

        assertTrue("No alignments read", count > 0);
    }

}
