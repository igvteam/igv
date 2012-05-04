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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam;

import org.broad.igv.sam.reader.GeraldParser;
import org.broad.tribble.readers.AsciiLineReader;
import org.junit.*;

import java.io.FileInputStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class GeraldParserTest {

    String geraldFile = "test/data/gerald/chr10_s_1_sorted.txt";

    public GeraldParserTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test parsing of Gerald record:
     * <p/>
     * HWUSI-EAS68R    0013    1       93      15585   5018    0       1       ACAGACTTCAAAAAACCGGAGCTTTTGCTGGGGATA
     * fffffffeffeeee^cccccfffffefeffffffff    chr10.fa                50497   R       36      34
     */
    @Test
    public void testParse() throws Exception {

        FileInputStream is = new FileInputStream(geraldFile);
        AsciiLineReader reader = new AsciiLineReader(is);
        GeraldParser instance = new GeraldParser();
        GeraldAlignment result = instance.readNextRecord(reader);

        String readName = "HWUSI-EAS68R:0013:1:93:15585:5018:0";
        assertEquals(readName, result.getReadName());
        assertTrue(result.isPaired());

        AlignmentBlock b = result.getAlignmentBlocks()[0];

        byte[] read = AlignmentUtils.reverseComplement("ACAGACTTCAAAAAACCGGAGCTTTTGCTGGGGATA").getBytes();
        for (int i = 0; i < read.length; i++) {
            assertEquals(read[i], b.getBases()[i]);
            byte q = b.getQualities()[i];
            assertTrue(q >= 0 && q <= 255);
        }

        assertEquals("chr10.fa", result.getChr());

        assertEquals(50496, result.getStart());
        assertEquals(50496 + read.length, result.getEnd());

        assertEquals(34, result.getMappingQuality());

    }

}