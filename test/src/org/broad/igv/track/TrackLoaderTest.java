package org.broad.igv.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static junit.framework.Assert.assertEquals;


/**
 * @author Jim Robinson
 * @date 10/3/11
 */
public class TrackLoaderTest {

    TrackLoader trackLoader;

    @Before
    public void setUp() throws Exception {
        TestUtils.setUpTestEnv();
        trackLoader = new TrackLoader();
    }

    @Test
    public void testLoadBEDIndexed() throws Exception {
        ResourceLocator locator = new ResourceLocator(TestUtils.DATA_DIR + "/bed/intervalTest.bed");
        if (!TrackLoader.isIndexed(locator.getPath())) {
            IgvTools igvTools = new IgvTools();
            igvTools.doIndex(locator.getPath(), IgvTools.LINEAR_INDEX, IgvTools.LINEAR_BIN_SIZE);
        }
        Genome genome = TestUtils.loadGenome();
        List<Track> tracks = trackLoader.load(locator, genome);
        assertEquals(1, tracks.size());
        Track track = tracks.get(0);
        assertEquals(locator, track.getResourceLocator());

    }

    @Test
    public void testLoadBEDNotIndexed() throws Exception {
        ResourceLocator locator = new ResourceLocator(TestUtils.DATA_DIR + "/bed/canFam2_alias.bed");
        if (TrackLoader.isIndexed(locator.getPath())) {
            File f = new File(locator.getPath() + ".idx");
            f.delete();
        }
        Genome genome = TestUtils.loadGenome();
        List<Track> tracks = trackLoader.load(locator, genome);
        assertEquals(1, tracks.size());
        Track track = tracks.get(0);
        assertEquals(locator, track.getResourceLocator());

    }

}
