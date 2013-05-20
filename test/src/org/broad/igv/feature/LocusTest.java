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

package org.broad.igv.feature;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class LocusTest {

    public LocusTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of isValid method, of class Locus.
     */
    @Test
    public void testValid() {
        String locusString = "chr1:15,000,000-20,000,000";
        Locus locus = Locus.fromString(locusString);
        assertNotNull(locus);
        assertTrue(locus.isValid());
        assertEquals("chr1", locus.getChr());
        assertEquals(15000000, locus.getStart());
        assertEquals(20000000, locus.getEnd());
    }

    public void testInvalid1() {
        String locusString = "somerandomestring";
        Locus locus = Locus.fromString(locusString);
        assertFalse(locus.isValid());
    }

    public void testInvalid2() {
        String locusString = "chr1:20-15";
        Locus locus = Locus.fromString(locusString);
        assertFalse(locus.isValid());
    }


}