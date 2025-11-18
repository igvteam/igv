package org.broad.igv.bedpe;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class HicSourceTest {

    static String testURL = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/cc0accc9-a89f-4b40-bac6-de754a7f77c9/4DNFILI6XOB5.hic";
    private static HicSource src;


    @BeforeClass
    public static void setUpClass() throws IOException {
        src = new HicSource(testURL,null);
    }


    @Test
    public void testGetFeatures() throws IOException {

        var features = src.getFeatures("1", 0, 248956422, 156872.35160680528, "NONE", 5000);
        assertEquals(50043, features.size());
    }


}