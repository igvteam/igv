/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.feature;

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
        String locusString = "chr1:15,000,001-20,000,000";
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