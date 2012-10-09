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

package org.broad.igv.dev.plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
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
 * User: jacob
 * Date: 2012-Oct-01
 */
public class AlignmentIOTest extends AbstractHeadlessTest {

    @Test
    public void testEncodeDecode() throws Exception {
        String testFile = TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam";
        SAMReader reader = new SAMReader(testFile);

        Iterator<Alignment> inputAlignmentIterator = reader.iterator();
        ArrayList<Alignment> inputAlignmentList = new ArrayList<Alignment>();

        while (inputAlignmentIterator.hasNext()) {
            Alignment al = inputAlignmentIterator.next();
            inputAlignmentList.add(al);
        }
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        FeatureEncoder<Alignment> alignmentEncoder = new SamAlignmentEncoder();
        alignmentEncoder.encodeAll(bos, reader.iterator());

        ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());

        FeatureDecoder<Alignment> alignmentDecoder = new AlignmentDecoder();
        Iterator<Alignment> decodedAlignments = alignmentDecoder.decodeAll(bis, false);

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
