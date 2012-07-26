package org.broad.igv.tdf;

import org.broad.igv.util.ResourceLocator;
import org.junit.Test;

import java.util.Set;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 7/19/12
 *         Time: 11:29 AM
 */
public class TDFReaderTest {

    @Test
    public void testReader() throws Exception {
        String url = "http://www.broadinstitute.org/igvdata/encode/hg18/broadHistone/SignalK562H3k4me3.tdf";

        TDFReader reader =  new TDFReader(new ResourceLocator(url));

        int version = reader.getVersion();
        assertEquals(3, version);
        assertTrue(reader.compressed);

        String [] trackNames = reader.getTrackNames();
        int nTracks = trackNames.length;
        assertEquals(1, nTracks);

        String trackName = trackNames[0];
        assertEquals(trackName, "SignalK562H3k4me3.wig.gz");

        Set<String> chrNames = reader.getChromosomeNames();
        int nChromosomes = chrNames.size();
        assertEquals(24, nChromosomes);

        String datasetName =  "/chr1/z0/mean";
        TDFDataset dataset = reader.getDataset(datasetName);
        assertEquals(datasetName, dataset.getName());

        TDFTile tile = reader.readTile(dataset, 0);
        assertNotNull(tile);

    }
}
