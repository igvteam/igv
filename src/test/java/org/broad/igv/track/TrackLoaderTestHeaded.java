package org.broad.igv.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * @author Jim Robinson
 * @date 10/3/11
 */
public class TrackLoaderTestHeaded extends AbstractHeadedTest {

    TrackLoader trackLoader;

    @Before
    public void setUp() throws Exception {
        super.setUp();
        trackLoader = new TrackLoader();
    }

    @Test
    public void testFilesHeaded() throws Exception {

        String ex_filename = "/vcf/example4-last-gsnap-2_fixed.vcf";
        Genome genome = TestUtils.loadGenome();
        List<String> finames = new ArrayList<String>(Arrays.asList(TrackLoaderTest.filenamesTryIndex));

        finames.add(ex_filename);

        for (String finame : finames) {
            TrackLoaderTest.tstLoadFi(trackLoader, TestUtils.DATA_DIR + finame, null, genome, true);
        }
    }


    @Test
    public void testReadVCF() throws Exception {
        String file = TestUtils.DATA_DIR + "vcf/outputPileup.flt1.vcf";
        TestUtils.createIndex(file);
        ResourceLocator locator = new ResourceLocator(file);
        //For files with 1 record, this threw a null pointer exception prior to r1595
        List<Track> tracks = igv.load(locator);

        Assert.assertEquals(1, tracks.size());
    }

}
