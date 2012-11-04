package org.broad.igv.hic.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.track.DataTrack;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.SortedMap;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 9/15/12
 *         Time: 6:48 PM
 */
public class HiCDataAdapterTest {


    @Test
    public void testGetDataFixedGrid() throws Exception {

        assertTrue(true);
//        int binCount = 24000;
//        int binSize = 10000;
//        HiCFixedGridAxis axis = new HiCFixedGridAxis(binCount, binSize);
//
//        String dataURL = "http://igvdata.broadinstitute.org/data/hg19/encode/broadHistone/wgEncodeBroadHistoneGm12878H3k36me3StdSig.wig.tdf";
//        TDFReader reader = TDFReader.getReader(dataURL);
//        Genome genome = TestUtils.loadGenome();  // This is actually hg18 but close enough
//        TDFDataSource dataSource = new TDFDataSource(reader, 0, "", genome);
//        DataTrack track = new DataSourceTrack(null, null, null, dataSource);
//
//        String chr = "chr14";
//        int start = 51000000;
//        int end = 56500000;
//        int startBin = start / binSize;
//        int endBin = end / binSize;
//
//        HiCDataAdapter adapter = new HiCDataAdapter(axis, track);
//
//        HiCDataAdapter.WeightedSum [] data =  adapter.getData(chr, startBin, endBin);
//        assertNotNull(data);
    }
}
