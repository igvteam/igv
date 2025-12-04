package org.broad.igv.feature.genome;

import static org.junit.Assert.*;

public class ClinVarTest {

    @org.junit.Test
    public void testClinVarRecordExists() {
        String existingHgvs = "NM_000546.5:c.215C>G"; // Known to exist
        String nonExistingHgvs = "NM_000000.0:c.9999A>T"; // Known not to exist
        assertNotNull(ClinVar.getClinVarURL(existingHgvs));
        assertNull(ClinVar.getClinVarURL(nonExistingHgvs));
    }
}