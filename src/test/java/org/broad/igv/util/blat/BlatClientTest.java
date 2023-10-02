package org.broad.igv.util.blat;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.Strand;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.util.MessageUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class BlatClientTest extends AbstractHeadlessTest  {

    @Test
    public void blatGET() throws IOException {

        //https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=GTCCTCGGAACCAGGACCTCGGCGTGGCCTAGCG&type=DNA&db=hg19&output=json
        //http://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=$SEQUENCE&type=DNA&db=$DB&output=json

        PreferencesManager.getPreferences().put(Constants.BLAT_URL, "https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=$SEQUENCE&type=DNA&db=$DB&output=json");
        List<PSLRecord> results = BlatClient.blat("hg19", "GTCCTCGGAACCAGGACCTCGGCGTGGCCTAGCG");

        assertEquals(1, results.size());
        assertEquals("chr21", results.get(0).getChr());
        assertEquals(Strand.POSITIVE, results.get(0).getStrand());
        assertEquals(34, results.get(0).getqSize());
    }

    @Test
    public void blatPOST() throws IOException {

        //https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=GTCCTCGGAACCAGGACCTCGGCGTGGCCTAGCG&type=DNA&db=hg19&output=json
        //http://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=$SEQUENCE&type=DNA&db=$DB&output=json

        PreferencesManager.getPreferences().put(Constants.BLAT_URL, "https://genome.ucsc.edu/cgi-bin/hgBlat");
        List<PSLRecord> results = BlatClient.blat("hg19", "GTCCTCGGAACCAGGACCTCGGCGTGGCCTAGCG");

        assertEquals(1, results.size());
        assertEquals("chr21", results.get(0).getChr());
        assertEquals(Strand.POSITIVE, results.get(0).getStrand());
        assertEquals(34, results.get(0).getqSize());
    }

    @Test
    public void blatTooShort() throws IOException {

        //https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=GTCCTCGGAACCAGGACCTCGGCGTGGCCTAGCG&type=DNA&db=hg19&output=json
        try {
            List<PSLRecord> results = BlatClient.blat("hg19", "GTCCTCGGA");
            fail("Exception expected");
        } catch (Exception e) {
            // This is expected
            assertTrue(true);

        }
    }


//    public static void main(String [] args) throws IOException {
//        (new BlatClientTest()).blatTooShort();
//    }


}