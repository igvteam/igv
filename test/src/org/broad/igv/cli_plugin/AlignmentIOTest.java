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

package org.broad.igv.cli_plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Iterator;

import static org.junit.Assert.assertEquals;

/**
 * Basic test to see if we can write out and read in alignments
 * @author jacob
 * @author 2012-Oct-01
 */
public class AlignmentIOTest extends AbstractHeadlessTest {

    @Test
    public void testEncodeDecode() throws Exception {
        String testFile = TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam";
        SAMReader reader = new SAMReader(testFile);

        Iterator<PicardAlignment> inputAlignmentIterator = reader.iterator();
        ArrayList<PicardAlignment> inputAlignmentList = new ArrayList<PicardAlignment>();

        while (inputAlignmentIterator.hasNext()) {
            PicardAlignment al = inputAlignmentIterator.next();
            inputAlignmentList.add(al);
        }
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        FeatureEncoder<PicardAlignment> alignmentEncoder = new SamAlignmentEncoder();
        alignmentEncoder.encodeAll(bos, reader.iterator());

        ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());

        FeatureDecoder<PicardAlignment> alignmentDecoder = new AlignmentDecoder();
        Iterator<PicardAlignment> decodedAlignments = alignmentDecoder.decodeAll(bis, false);

        int ind = 0;
        while (decodedAlignments.hasNext()) {

            Alignment act = decodedAlignments.next();
            Alignment exp = inputAlignmentList.get(ind++);
            TestUtils.assertFeaturesEqual(exp, act);
            assertEquals(exp.getCigarString(), act.getCigarString());
        }

        assertEquals("Different number of alignments read in as out", inputAlignmentList.size(), ind);
    }
}
