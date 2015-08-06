/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.sam.reader;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.tools.IGVToolsTest;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * @author jacob
 * @since 2012/01/25
 */
public class MergedAlignmentReaderTest extends AbstractHeadlessTest {

    public static String[] generateRepLargebamsList(File listFile) throws IOException{
        listFile.delete();
        listFile.deleteOnExit();
        String listPath = listFile.getPath();

        return IGVToolsTest.generateRepLargebamsList(listPath, "HG00171.hg18.bam", 2);
    }

    @Test
    public void testSimpleRead() throws Exception {

        File listFile = new File(TestUtils.LARGE_DATA_DIR, "2largebams.bam.list");

        String[] actfiles = generateRepLargebamsList(listFile);
        int start = 151667156;
        int end = start + 10000;
        int num_combined = 0;
        int num_sep = 0;
        Alignment align;

        AlignmentReader mergedReader = AlignmentReaderFactory.getBamListReader(listFile.getAbsolutePath(), false);
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

        BufferedReader in = new BufferedReader(new FileReader(listFile));
        String singfile = in.readLine();
        singfile = FileUtils.getAbsolutePath(singfile, listFile.getPath());
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
