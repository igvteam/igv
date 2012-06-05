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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.tools.IGVToolsTest;
import org.broad.igv.util.TestUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.util.Iterator;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012/03/01
 */
public class AlignmentReaderFactoryTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 1e3 * 60);

    @Test
    public void testGetMergedReader() throws Exception {
        String relfiname = "HG00171.hg18.bam";
        String path = TestUtils.LARGE_DATA_DIR + "temprel.bam.list";
        String[] relnames = IGVToolsTest.generateRepLargebamsList(path, relfiname, 3, false);
        for (String s : relnames) {
            assertFalse((new File(s)).isAbsolute());
        }

        String apath = TestUtils.LARGE_DATA_DIR + "tempabs.bam.list";
        String[] absnames = IGVToolsTest.generateRepLargebamsList(apath, relfiname, 3, true);
        for (String s : absnames) {
            assertTrue((new File(s)).isAbsolute());
        }

        AlignmentReader relreader = AlignmentReaderFactory.getReader(path, false);
        AlignmentReader absreader = AlignmentReaderFactory.getReader(apath, false);
        Iterator<Alignment> iter = relreader.iterator();
        Iterator<Alignment> absiter = absreader.iterator();

        while (iter.hasNext()) {
            Alignment relA = iter.next();
            Alignment absA = absiter.next();
            //Only check 1 field for speed, these are long strings so unlikely to be equal by accident
            assertEquals(absA.getCigarString(), relA.getCigarString());
        }
    }

    public static void main(String[] args) {
        String path = "http://www.example.com";
        File fi = new File(path);
        System.out.println(fi.isAbsolute());
    }


}
