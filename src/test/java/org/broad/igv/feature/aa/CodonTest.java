package org.broad.igv.feature.aa;

import htsjdk.tribble.NamedFeature;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.IGVNamedFeature;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012-Oct-04
 */
public class CodonTest extends AbstractHeadlessTest {

    private BasicFeature egfr;

    @Override
    public void setUp() throws Exception {
        super.setUp();
        egfr = (BasicFeature) genome.getFeatureDB().getFeature("egfr");
    }

    @Test @Ignore("Fails unless tests are run in separate JVMs")
    public void testGetCodon() {
        // See http://www.ncbi.nlm.nih.gov/nuccore/NM_201283
        //Note: This covers a break in exons
        char[] expected = "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLE".toCharArray();
        tstGetCodon(egfr, expected);

    }

    @Test @Ignore("Fails unless tests are run in separate JVMs")
    public void testGetCodonNeg() {
        //See http://www.ncbi.nlm.nih.gov/nuccore/NM_004985
        String exp_string = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV";
        char[] expected = exp_string.toCharArray();

        BasicFeature KRAS = (BasicFeature) genome.getFeatureDB().getFeature("KRAS");
        tstGetCodon(KRAS, expected);
    }

    public void tstGetCodon(IGVNamedFeature feature, char[] expected) {
        BasicFeature bf = (BasicFeature) feature;

        for (int pos = 0; pos < expected.length; pos++) {
            Codon codon = bf.getCodon(genome, bf.getChr(), pos + 1);
            assertEquals(expected[pos], codon.getAminoAcid().getSymbol());
        }
    }

    //Check that there is no exception when loading an
    @Test
    public void testGetCodonInvalid() {

        String geneName = "APOL1";
        int proteinPos = 384;
        //The position isn't valid for this transcript only
        String invalidId = "NM_001136541";
        //valid for these ids:
        //NM_003661, NM_145343, NM_001136540

        List<NamedFeature> featuresList = genome.getFeatureDB().getFeaturesMatching(geneName);
        for (NamedFeature feat : featuresList) {
            BasicFeature bf = (BasicFeature) feat;
            System.out.println(bf.getIdentifier());
            Codon c = bf.getCodon(genome, bf.getChr(),  proteinPos);

            if (bf.getIdentifier().equals(invalidId)) {
                assertNull(c);
            } else {
                assertNotNull(c);
                assertTrue(c.isGenomePositionsSet());
            }


        }

    }
}
