/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.tools.IGVToolsTest;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
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
public class MergedAlignmentReaderTest {


    @Before
    public void setUp() throws Exception {
        Globals.setHeadless(true);

    }

    @After
    public void tearDown() throws Exception {

    }

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
        SAMFileHeader mergedHeader = mergedReader.getHeader();
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
        SAMFileHeader singHeader = singReader.getHeader();

        //Compare the headers
        assertEquals(singHeader.getTextHeader(), mergedHeader.getTextHeader());
        assertEquals(singHeader.getVersion(), mergedHeader.getVersion());
        assertEquals(singHeader.getCreator(), mergedHeader.getCreator());


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
