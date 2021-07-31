package org.broad.igv.feature.aa;

import org.junit.Ignore;
import org.junit.Test;
import org.junit.Before;

import static org.junit.Assert.*;

public class CodonTableManagerTest {

    @Before
    public void setup() {
        CodonTableManager.getInstance().resetToDefaults();
    }

    @Test
    public void getCodonTable() {

        CodonTable table = CodonTableManager.getInstance().getCodonTableForChromosome("hg19", "chr1");
        assertEquals(1, table.getId());

        table = CodonTableManager.getInstance().getCodonTableForChromosome("hg19", "chrM");
        assertEquals(2, table.getId());
    }


    //Just check that we have the ids we expect
    @Test
    public void testCodonTablesExist() throws Exception {
        int[] expIds = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24};
        for (int id : expIds) {
            assertNotNull(CodonTableManager.getInstance().getCodonTableByID(id));
        }
    }


    //Sample a few codon tables, test different translations
    @Test
    public void testVariousCodonTables()  {

        int[] codonTableIds = {2, 2, 3, 6, 16, 22, 23, 24};
        String[] testCodons = {"AGA", "TGA", "CTG", "TAA", "TAG", "TCA", "TTA", "AGA"};
        char[] expAAs = {'*', 'W', 'T', 'Q', 'L', '*', '*', 'S'};

        for (int ii = 0; ii < expAAs.length; ii++) {
            int id = codonTableIds[ii];
            CodonTable codonTable = CodonTableManager.getInstance().getCodonTableByID(id);
            AminoAcid actualAA = codonTable.getAminoAcid(testCodons[ii]);
            assertNotSame("Got null amino acid for " + testCodons[ii], AminoAcid.NULL_AMINO_ACID, actualAA);
            assertEquals(String.valueOf(expAAs[ii]), String.valueOf(actualAA.getSymbol()));

            //We want to only store one copy of each amino acid. Check that this is the case
            assertTrue(AminoAcidManager.getAminoAcidByName(actualAA.getShortName()) == actualAA);
        }
    }

//    public static void main(String [] args) {
//        (new CodonTableManagerTest()).testVariousCodonTables();
//    }
}