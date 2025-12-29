package org.igv.feature.genome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

public class ClinVarTest {

    @org.junit.Test
    public void testClinVarRecordExists() {
        String existingHgvs = "NM_000546.5:c.215C>G"; // Known to exist
        String nonExistingHgvs = "NM_000000.0:c.9999A>T"; // Known not to exist
        String clinvarURL = ClinVar.getClinVarURL(existingHgvs);
        String expectedValue = "https://www.ncbi.nlm.nih.gov/clinvar/variation/237944/";
        assertEquals(expectedValue, clinvarURL);

        assertNull(ClinVar.getClinVarURL(nonExistingHgvs));
    }
}