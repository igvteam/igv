package org.broad.igv.ga4gh;

import com.google.gson.JsonObject;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class Ga4ghAPIHelperTest {

    static Ga4ghProvider provider = Ga4ghAPIHelper.GA4GH_GOOGLE_PROVIDER;


    @Test
    public void testReadGroupsetsSearch() throws Exception {

        String datasetId = "10473108253681171589";

        List<Ga4ghReadset> readGroupsets = Ga4ghAPIHelper.searchReadGroupsets(provider, datasetId, 1000);

        assertNotNull(readGroupsets);

    }

    @Test
    public void testReferenceSearch() throws Exception {

        String referenceSetId = "EOSt9JOVhp3jkwE";

        List<JsonObject> references = Ga4ghAPIHelper.searchReferences(provider, referenceSetId, 1000);

        assertNotNull(references);

    }


    @Test
    public void testReads() throws Exception {

        String readGroupsetId = "CMvnhpKTFhCwvIWYw9eikzQ";
    }
}