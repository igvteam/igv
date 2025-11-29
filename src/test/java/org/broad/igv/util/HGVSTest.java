package org.broad.igv.util;

import org.junit.Test;

import static org.junit.Assert.*;

public class HGVSTest {

    @Test
    public void testIsValidHGVS() {
        //assertTrue(HGVS.isValidHGVS("NM_000546.5:c.215C>G"));
        //assertTrue(HGVS.isValidHGVS("ENST00000380152.6:c.215C>G"));
        assertTrue(HGVS.isValidHGVS("NC_000017.11:g.7579472C>G"));
        assertTrue(HGVS.isValidHGVS("NC_000017.11:g.7579472"));
        assertFalse(HGVS.isValidHGVS("Invalid_HGVS_String"));
        assertFalse(HGVS.isValidHGVS("NC_000017.11:g.7579472C>"));
    }
}