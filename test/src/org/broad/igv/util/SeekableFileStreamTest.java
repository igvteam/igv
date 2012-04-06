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

package org.broad.igv.util;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.SeekableFileStream;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 6, 2009
 * Time: 1:07:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class SeekableFileStreamTest {


    @Test
    public void testSeek() throws Exception {
        String expectedLine = "chr22\t14000000\t15000000\trecombRate\t0.182272\t0.20444\t0.160104\t0\t0\t0\t0\t0\t0";
        File testFile = new File(TestUtils.DATA_DIR + "igv/recombRate.igv.txt");
        SeekableFileStream is = new SeekableFileStream(testFile);
        is.seek(149247);
        AsciiLineReader reader = new AsciiLineReader(is);
        String nextLine = reader.readLine();
        assertEquals(expectedLine, nextLine);
    }
}
