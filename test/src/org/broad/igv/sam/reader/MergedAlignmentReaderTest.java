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

package org.broad.igv.sam.reader;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.tools.IGVToolsTest;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/01/25
 */
public class MergedAlignmentReaderTest extends AbstractHeadlessTest {

    @Test
    public void testSimpleRead() throws Exception {
        //This test file is 2 lines of the same file.
        //Check that we get the same results twice
        String file = TestUtils.LARGE_DATA_DIR + File.separator + "2largebams.bam.list";
        String[] actfiles = IGVToolsTest.generateRepLargebamsList(file, "HG00171.hg18.bam", 2);
        int start = 151667156;
        int end = start + 10000;
        int num_combined = 0;
        int num_sep = 0;
        Alignment align;

        AlignmentReader mergedReader = AlignmentReaderFactory.getBamListReader(file, false);
        CloseableIterator<Alignment> combData = mergedReader.query("chr1", start, end, false);

        Map<Float, Integer> combinedCounts = new HashMap();
        while (combData.hasNext()) {
            align = combData.next();
            num_combined += 1;
            float score = align.getScore();
            Integer oldcount = combinedCounts.get(score);
            int newcount = oldcount == null ? 1 : oldcount + 1;
            combinedCounts.put(score, newcount);
        }
        mergedReader.close();

        BufferedReader in = new BufferedReader(new FileReader(file));
        String singfile = in.readLine();
        singfile = FileUtils.getAbsolutePath(singfile, file);
        AlignmentReader singReader = AlignmentReaderFactory.getReader(singfile, false);

        CloseableIterator<Alignment> singData = singReader.query("chr1", start, end, false);

        Map<Float, Integer> singCounts = new HashMap();
        while (singData.hasNext()) {
            align = singData.next();
            num_sep += 1;
            float score = align.getScore();
            Integer oldcount = singCounts.get(score);
            int newcount = oldcount == null ? 1 : oldcount + 1;
            singCounts.put(score, newcount);
        }

        singReader.close();

        assertEquals(2 * num_sep, num_combined);
        assertEquals(singCounts.size(), combinedCounts.size());
        for (Float score : singCounts.keySet()) {
            assertEquals(2 * singCounts.get(score), (int) combinedCounts.get(score));
        }
    }

}
